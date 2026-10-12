#!/bin/bash
set -euo pipefail

LIU_DIR="$HOME/sud/liu-2019"
OUT_DIR="$HOME/sud/liu-2019/liftover_only"
CHAIN="$HOME/liftover_files/hg19ToHg38.over.chain.gz"
LIFTOVER_BIN="$HOME/liftOver"

mkdir -p "$OUT_DIR"

if [ ! -x "$LIFTOVER_BIN" ]; then
    echo "ERROR: liftOver binary not found or not executable: $LIFTOVER_BIN"
    echo "Try: chmod +x $LIFTOVER_BIN"
    exit 1
fi

if [ ! -f "$CHAIN" ]; then
    echo "ERROR: chain file not found: $CHAIN"
    exit 1
fi

cat > /tmp/liftover_liu_only.py << 'EOF'
import subprocess
import pandas as pd
from pathlib import Path

liu_dir = Path("~/sud/liu-2019").expanduser()
out_dir = Path("~/sud/liu-2019/liftover_only").expanduser()
chain = Path("~/liftover_files/hg19ToHg38.over.chain.gz").expanduser()
liftover_bin = Path("~/liftOver").expanduser()

out_dir.mkdir(parents=True, exist_ok=True)

files = [
    "AgeOfInitiation.WithoutUKB.txt",
    "AgeofInitiation.txt",
    "CigarettesPerDay.WithoutUKB.txt",
    "CigarettesPerDay.txt",
    "DrinksPerWeek.WithoutUKB.txt",
    "DrinksPerWeek.txt",
    "SmokingCessation.WithoutUKB.txt",
    "SmokingCessation.txt",
    "SmokingInitiation.WithoutUKB.txt",
    "SmokingInitiation.txt"
]

def is_snv(a):
    return isinstance(a, str) and len(a) == 1 and a.upper() in {"A", "C", "G", "T"}

def is_ambiguous(a1, a2):
    pair = {str(a1).upper(), str(a2).upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}

for filename in files:
    print("Processing", filename)

    in_file = liu_dir / filename

    df = pd.read_csv(
        in_file,
        sep=r"\s+",
        dtype={"RSID": str},
        low_memory=False
    )

    df.columns = df.columns.str.strip()

    df["REF"] = df["REF"].str.upper()
    df["ALT"] = df["ALT"].str.upper()

    df = df[
        df["REF"].apply(is_snv) &
        df["ALT"].apply(is_snv)
    ].copy()

    df = df[
        ~df.apply(lambda r: is_ambiguous(r["REF"], r["ALT"]), axis=1)
    ].copy()

    df["CHROM"] = df["CHROM"].astype(str)
    df["POS"] = df["POS"].astype(int)

    df["row_id"] = ["var{}".format(i) for i in range(len(df))]

    pre_file = out_dir / "{}.pre_liftover.tsv".format(filename)
    bed_hg37 = out_dir / "{}.hg37.bed".format(filename)
    bed_hg38 = out_dir / "{}.hg38.bed".format(filename)
    unmapped = out_dir / "{}.unmapped.bed".format(filename)
    lifted_tsv = out_dir / "{}.lifted_hg38.tsv".format(filename)

    df.to_csv(pre_file, sep="\t", index=False)

    bed = pd.DataFrame({
        "chrom": "chr" + df["CHROM"],
        "start": df["POS"] - 1,
        "end": df["POS"],
        "name": df["row_id"]
    })

    bed.to_csv(
        bed_hg37,
        sep="\t",
        header=False,
        index=False
    )

    subprocess.run(
        [
            str(liftover_bin),
            str(bed_hg37),
            str(chain),
            str(bed_hg38),
            str(unmapped)
        ],
        check=True
    )

    lifted = pd.read_csv(
        bed_hg38,
        sep="\t",
        header=None,
        names=["chr_hg38", "start_hg38", "end_hg38", "row_id"],
        dtype={"row_id": str}
    )

    lifted["CHR_hg38"] = lifted["chr_hg38"].str.replace("chr", "", regex=False)
    lifted["BP_hg38"] = lifted["end_hg38"].astype(int)

    out = df.merge(
        lifted[["row_id", "CHR_hg38", "BP_hg38"]],
        on="row_id",
        how="inner"
    )

    out.to_csv(
        lifted_tsv,
        sep="\t",
        index=False
    )

    print("  Variants prepared: {:,}".format(len(df)))
    print("  Lifted variants: {:,}".format(len(out)))
    print("  Wrote:", lifted_tsv)
    print("  Unmapped:", unmapped)

print("Done.")
EOF

python3 /tmp/liftover_liu_only.py
