#!/usr/bin/env bash
# For each site:
#   - Slim log:  logs/<SITE>_<timestamp>.log	w/ some of the more important info
#   - Full log:  logs/<SITE>_<timestamp>.full.log.gz	w/ the full log but compressed
# Appends to:    /home/nvalenci/polmm/maf1e-2_site_completion.txt
# Adds column to *completion.txt file:   err_non_numeric_matrix_extent = Yes/No if error text is found in full log.

set -uo pipefail

# Fixed paths
WORKDIR="/home/nvalenci/polmm"
OUTPUT_DIR="${WORKDIR}/output/Nas0_MAF1e-2"
LOGDIR="${WORKDIR}/logs"
COMPLETION_FILE="${WORKDIR}/maf1e-2_site_completion.txt"

# declaring which sites and .sh files to run
declare -a SITES=("AAU" "KEMRI" "Uganda" "Moi" "UCT")
declare -a SCRIPTS=("aau_run_polmm.sh" "kemri_run_polmm.sh" "uganda_run_polmm.sh" "moi_run_polmm.sh" "uct_run_polmm.sh")

mkdir -p "${OUTPUT_DIR}" "${LOGDIR}"

# Create header if missing (includes error column)
if [[ ! -f "${COMPLETION_FILE}" ]]; then
  printf "timestamp\tsite\tscript_exit_code\toutputs_found\tmatching_files\terr_non_numeric_matrix_extent\n" \
    > "${COMPLETION_FILE}"
fi

# Detect expected outputs (supports both flat and subdir layouts)
check_outputs() {
  local site="$1"
  local found="no"
  local matched_list="-"

  local p1="${OUTPUT_DIR}/assoc_${site}_Nas0_MAF1e-2"* 
  local p2="${OUTPUT_DIR}/Nas0_MAF1e-2/assoc_${site}_Nas0_MAF1e-2"*

  shopt -s nullglob
  local matches=()
  for p in "$p1" "$p2"; do
    for f in $p; do
      # only non-empty regular files
      if [[ -f "$f" && -s "$f" ]]; then
        matches+=("$f")
      fi
    done
  done
  shopt -u nullglob

  if (( ${#matches[@]} > 0 )); then
    found="yes"
    local rels=()
    for f in "${matches[@]}"; do rels+=("${f#${WORKDIR}/}"); done
    matched_list=$(IFS=, ; echo "${rels[*]}")
  fi

  echo "$found|$matched_list"
}

run_one() {
  local site="$1"
  local script="$2"

  local ts rc outputs found list errflag
  ts="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "==== [${ts}] Starting ${site} via ${script} ===="

  chmod +x "${WORKDIR}/${script}" 2>/dev/null || true

  # Per-site logs
  local runid
  runid="$(date -u +%Y%m%dT%H%M%SZ)"
  local slimlog="${LOGDIR}/${site}_${runid}.log"
  local fulllog="${LOGDIR}/${site}_${runid}.full.log.gz"
  local tmppath
  tmppath="$(mktemp)"

  # Run and capture ALL stdout+stderr to a temp file
  (
    cd "${WORKDIR}" && "./${script}"
  ) > "${tmppath}" 2>&1
  rc=$?

  # Build a compact, human-readable log:
  # Keep our headings, completion lines, errors, warnings, and any "Output:" echoes.
  # (Adjust the grep to include/exclude more patterns as you like.)
  grep -Ei \
    '^(====|\-\-\-\-|Error|Warning|Analysis Complete|POLMM association test completed|Output:)' \
    "${tmppath}" > "${slimlog}" || true

  # Save the full, unfiltered log compressed
  gzip -9c "${tmppath}" > "${fulllog}"
  rm -f "${tmppath}"

  # This 'Error in matrix..." error has been the very same thing we have been erroring with
  # for specific sites
	# trying to log which sites get this error
  if zgrep -q 'Error in matrix(adjGVec, n1, J - 1) : non-numeric matrix extent' "${fulllog}"; then
    errflag="Yes"
  else
    errflag="No"
  fi

  # Detect outputs
  outputs="$(check_outputs "${site}")"
  found="${outputs%%|*}"
  list="${outputs#*|}"

  # Append one TSV line
  printf "%s\t%s\t%d\t%s\t%s\t%s\n" "${ts}" "${site}" "${rc}" "${found}" "${list}" "${errflag}" \
    >> "${COMPLETION_FILE}"

  echo "---- ${site}: exit=${rc}, outputs_found=${found}, err_non_numeric_matrix_extent=${errflag}"
  echo "     Slim log: ${slimlog}"
  echo "     Full log: ${fulllog}"
  [[ "${found}" == "yes" ]] && echo "     Files: ${list}"
  echo
}

# Main sequence
for i in "${!SITES[@]}"; do
  run_one "${SITES[$i]}" "${SCRIPTS[$i]}"
done

echo "Done. Completion log: ${COMPLETION_FILE}"
