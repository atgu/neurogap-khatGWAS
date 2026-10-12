#!/usr/bin/env bash
set -euo pipefail

JOHNSON_DIR="${HOME}/sud/johnson-2020"
OUT_DIR="${HOME}/sud/johnson-2020/liftover_only"
CHAIN="${HOME}/liftover_files/hg19ToHg38.over.chain.gz"
LIFTOVER_BIN="${HOME}/liftOver"

mkdir -p "${OUT_DIR}"

[[ -x "${LIFTOVER_BIN}" ]] || { echo "ERROR: liftOver not executable: ${LIFTOVER_BIN}" >&2; exit 1; }
[[ -f "${CHAIN}" ]] || { echo "ERROR: chain file not found: ${CHAIN}" >&2; exit 1; }

cat > /tmp/liftover_johnson_effect_a1.py <<'PYTHON'
import subprocess
import pandas as pd
from pathlib import Path

johnson_dir = Path("~/sud/johnson-2020").expanduser()
out_dir = Path("~/sud/johnson-2020/liftover_only").expanduser()
chain = Path("~/liftover_files/hg19ToHg38.over.chain.gz").expanduser()
liftover_bin = Path("~/liftOver").expanduser()

files = [
    "CUD_AFR_full_public_11.14.2020",
    "CUD_EUR_casecontrol_public_11.14.2020",
    "CUD_EUR_full_public_11.14.2020",
]

def is_snv(a):
    return isinstance(a, str) and len(a) == 1 and a.upper() in {"A", "C", "G", "T"}

def is_ambiguous(a1, a2):
    return {a1, a2} in ({"A", "T"}, {"C", "G"})

for filename in files:
    print(f"Processing {filename}", flush=True)
    in_file = johnson_dir / filename
    if not in_file.exists():
        raise FileNotFoundError(in_file)

    df = pd.read_csv(in_file, sep=r"\s+", dtype={"SNP": str}, low_memory=False)
    df.columns = df.columns.str.strip()

    required = {"CHR", "SNP", "BP", "A1", "A2"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{in_file}: missing columns: {', '.join(missing)}")
    if "Beta" not in df.columns and "Z" not in df.columns:
        raise ValueError(f"{in_file}: neither Beta nor Z is present")

    # Johnson documentation: A1 is the effect allele and A2 is non-effect.
    # Liftover changes coordinates only; A1, A2, Beta/Z, and all other fields
    # remain attached to the same original row.
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
    pd.DataFrame({
        "chrom": "chr" + df["CHR"],
        "start": df["BP"] - 1,
        "end": df["BP"],
        "name": df["row_id"],
    }).to_csv(bed_hg19, sep="\t", header=False, index=False)

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

    out = df.merge(lifted[["row_id", "CHR_hg38", "BP_hg38"]], on="row_id", how="inner")
    out.to_csv(lifted_tsv, sep="\t", index=False)

    print(f"  Filtered variants: {len(df):,}")
    print(f"  Lifted variants:   {len(out):,}")
    print(f"  Wrote: {lifted_tsv}")

print("Done.", flush=True)
PYTHON

python3 /tmp/liftover_johnson_effect_a1.py
