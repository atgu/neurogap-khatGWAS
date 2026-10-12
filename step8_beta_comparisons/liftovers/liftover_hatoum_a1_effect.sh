#!/usr/bin/env bash
set -euo pipefail

HATOUM_DIR="${HOME}/sud/hatoum-2023"
OUT_DIR="${HOME}/sud/hatoum-2023/liftover_only"
CHAIN="${HOME}/liftover_files/hg19ToHg38.over.chain.gz"
LIFTOVER_BIN="${HOME}/liftOver"

mkdir -p "${OUT_DIR}"

[[ -x "${LIFTOVER_BIN}" ]] || {
  echo "ERROR: liftOver binary not found or not executable: ${LIFTOVER_BIN}" >&2
  exit 1
}
[[ -f "${CHAIN}" ]] || {
  echo "ERROR: chain file not found: ${CHAIN}" >&2
  exit 1
}

cat > /tmp/liftover_hatoum_effect_a1.py <<'PYTHON'
import subprocess
import pandas as pd
from pathlib import Path

hatoum_dir = Path("~/sud/hatoum-2023").expanduser()
out_dir = Path("~/sud/hatoum-2023/liftover_only").expanduser()
chain = Path("~/liftover_files/hg19ToHg38.over.chain.gz").expanduser()
liftover_bin = Path("~/liftOver").expanduser()

out_dir.mkdir(parents=True, exist_ok=True)

files = [
    "Hatoum2023AddictionAfrican.txt",
    "Hatoum2023AddictionEuropean.txt",
]

def is_snv(a):
    return isinstance(a, str) and len(a) == 1 and a.upper() in {"A", "C", "G", "T"}

def is_ambiguous(a1, a2):
    return {a1, a2} in ({"A", "T"}, {"C", "G"})

for filename in files:
    print(f"Processing {filename}", flush=True)
    in_file = hatoum_dir / filename
    if not in_file.exists():
        raise FileNotFoundError(in_file)

    df = pd.read_csv(in_file, sep=r"\s+", dtype={"SNP": str}, low_memory=False)
    df.columns = df.columns.str.strip()

    if "CHR" not in df.columns and "Chr" in df.columns:
        df = df.rename(columns={"Chr": "CHR"})

    required = {"CHR", "SNP", "BP", "A1", "A2", "P"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{in_file}: missing required columns: {', '.join(missing)}")

    if "European" in filename and "Beta" not in df.columns:
        raise ValueError(f"{in_file}: expected European Beta column")
    if "African" in filename and "OR" not in df.columns:
        raise ValueError(f"{in_file}: expected African OR column")

    # Hatoum documentation:
    #   A1 = effect allele; A2 = other allele.
    #   European effect = Beta; African effect = OR.
    # Liftover changes coordinates only. It does not swap alleles or modify
    # Beta/OR, so every effect remains attached to its documented A1 allele.
    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()

    df = df[
        df["SNP"].notna()
        & ~df["SNP"].isin([".", "NA"])
        & df["A1"].apply(is_snv)
        & df["A2"].apply(is_snv)
    ].copy()
    df = df[~df.apply(lambda r: is_ambiguous(r["A1"], r["A2"]), axis=1)].copy()

    df["CHR"] = df["CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    df["BP"] = pd.to_numeric(df["BP"], errors="raise").astype(int)
    df["row_id"] = [f"var{i}" for i in range(len(df))]

    pre_file = out_dir / f"{filename}.pre_liftover.tsv"
    bed_hg19 = out_dir / f"{filename}.hg19.bed"
    bed_hg38 = out_dir / f"{filename}.hg38.bed"
    unmapped = out_dir / f"{filename}.unmapped.bed"
    lifted_tsv = out_dir / f"{filename}.lifted_hg38.tsv"

    df.to_csv(pre_file, sep="\t", index=False)

    bed = pd.DataFrame({
        "chrom": "chr" + df["CHR"],
        "start": df["BP"] - 1,
        "end": df["BP"],
        "name": df["row_id"],
    })
    bed.to_csv(bed_hg19, sep="\t", header=False, index=False)

    subprocess.run(
        [str(liftover_bin), str(bed_hg19), str(chain), str(bed_hg38), str(unmapped)],
        check=True,
    )

    lifted = pd.read_csv(
        bed_hg38,
        sep="\t",
        header=None,
        names=["chr_hg38", "start_hg38", "end_hg38", "row_id"],
        dtype={"row_id": str},
    )
    lifted["CHR_hg38"] = lifted["chr_hg38"].str.replace("chr", "", regex=False)
    lifted["BP_hg38"] = lifted["end_hg38"].astype(int)

    out = df.merge(
        lifted[["row_id", "CHR_hg38", "BP_hg38"]],
        on="row_id",
        how="inner",
    )
    out.to_csv(lifted_tsv, sep="\t", index=False)

    print(f"  Input variants after filtering: {len(df):,}")
    print(f"  Lifted variants:               {len(out):,}")
    print(f"  Wrote: {lifted_tsv}")
    print(f"  Unmapped: {unmapped}")

print("Done.", flush=True)
PYTHON

python3 /tmp/liftover_hatoum_effect_a1.py
