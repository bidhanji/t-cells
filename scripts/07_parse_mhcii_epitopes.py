"""
Parse DTU NetMHCIIpan-4.3 batch output files.

IMPORTANT: These files have a .xls extension but are NOT real Excel binary
files and are NOT HTML. They are plain tab-separated text with a two-row
header (allele names on row 2, sub-column names on row 3). Do not feed them
to pd.read_excel() or pd.read_html() -- they will fail or silently mis-parse.

File layout (0-indexed lines):
  line 0: the NetMHCIIpan command that was run (metadata, discarded)
  line 1: allele names, each appearing once above its 4-column block
  line 2: sub-column headers: Pos, Peptide, ID, Target,
          then (Core, Inverted, Score_EL, Rank_EL) x N_alleles,
          then Ave, NB
  line 3+: tab-separated data rows
"""

import glob
import os
import sys

import pandas as pd

# ---- config -----------------------------------------------------------
# INPUT_DIR is resolved relative to THIS SCRIPT'S location, not whatever
# folder you happened to launch python from -- that mismatch is exactly
# what caused "no files matching" last time. Set this to wherever your
# .xls files actually sit relative to the scripts/ folder.
# Examples: ".." (project root), "../data", "../data/mhcii_raw", "."
INPUT_DIR = "../data/processed/batches"
FILE_PATTERN = "dtu_mhcii_batch_*.xls"
OUTPUT_CSV = "mhcii_strong_binders.csv"
RANK_STRONG_THRESHOLD = 1.0         # must match -rankS used at submission
MIN_ALLELES_HIT = 1                 # NB >= this many alleles counts as "strong binder" row
FIXED_COLS_EXPECTED = ["Pos", "Peptide", "ID", "Target"]
METRIC_BLOCK = ["Core", "Inverted", "Score_EL", "Rank_EL"]
TAIL_COLS_EXPECTED = ["Ave", "NB"]
# ------------------------------------------------------------------------


def parse_dtu_mhcii_file(filepath: str) -> pd.DataFrame:
    """Parse one DTU NetMHCIIpan .xls (really TSV) file into a flat DataFrame."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    if len(lines) < 4:
        raise ValueError(
            f"File has only {len(lines)} lines -- expected a command line, "
            f"allele line, sub-header line, and >=1 data row. This file is "
            f"probably empty, truncated, or not a real NetMHCIIpan export."
        )

    command_line = lines[0].strip()
    if "NetMHCIIpan" not in command_line and "netMHCIIpan" not in command_line:
        raise ValueError(
            f"Line 1 does not look like a NetMHCIIpan command line: "
            f"{command_line[:80]!r}. Refusing to guess the layout -- "
            f"open this file manually and check what it actually contains."
        )

    allele_line = lines[1].rstrip("\n").split("\t")
    alleles = [a.strip() for a in allele_line if a.strip()]
    if len(alleles) == 0:
        raise ValueError("Line 2 (allele row) is empty -- can't determine allele panel.")

    subheader = lines[2].rstrip("\n").split("\t")

    fixed_cols = subheader[:4]
    if fixed_cols != FIXED_COLS_EXPECTED:
        raise ValueError(
            f"Sub-header's first 4 columns are {fixed_cols}, expected "
            f"{FIXED_COLS_EXPECTED}. File format has changed or this isn't "
            f"a NetMHCIIpan -xls export -- do not proceed blind."
        )

    tail_cols = subheader[-2:]
    if tail_cols != TAIL_COLS_EXPECTED:
        raise ValueError(
            f"Sub-header's last 2 columns are {tail_cols}, expected "
            f"{TAIL_COLS_EXPECTED}. Aborting rather than mis-assigning columns."
        )

    metric_block = subheader[4:-2]
    expected_metric_len = len(alleles) * len(METRIC_BLOCK)
    if len(metric_block) != expected_metric_len:
        raise ValueError(
            f"Found {len(alleles)} alleles on line 2 but {len(metric_block)} "
            f"metric columns on line 3 (expected {expected_metric_len} = "
            f"{len(alleles)} alleles x {len(METRIC_BLOCK)} metrics each). "
            f"Allele count and column count don't agree -- do not proceed."
        )

    # sanity check each 4-column block actually matches Core/Inverted/Score_EL/Rank_EL
    for i in range(len(alleles)):
        block = metric_block[i * 4:(i + 1) * 4]
        if block != METRIC_BLOCK:
            raise ValueError(
                f"Metric block for allele #{i+1} ({alleles[i]}) is {block}, "
                f"expected {METRIC_BLOCK}. Column order assumption is wrong."
            )

    cols = list(fixed_cols)
    for allele in alleles:
        cols += [f"{allele}_{m}" for m in METRIC_BLOCK]
    cols += tail_cols

    df = pd.read_csv(
        filepath,
        sep="\t",
        skiprows=3,
        header=None,
        names=cols,
        engine="python",
    )

    if df.empty:
        raise ValueError("Parsed 0 data rows -- file has header only, no peptides.")

    df["_source_file"] = os.path.basename(filepath)
    df["_alleles"] = ",".join(alleles)
    return df


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_dir_resolved = os.path.abspath(os.path.join(script_dir, INPUT_DIR))
    search_glob = os.path.join(input_dir_resolved, FILE_PATTERN)
    filepaths = sorted(glob.glob(search_glob))

    if not filepaths:
        print(f"ERROR: no files matching {FILE_PATTERN!r} found.")
        print(f"  Script location:        {script_dir}")
        print(f"  Resolved search folder: {input_dir_resolved}")
        print(f"  (INPUT_DIR is currently set to {INPUT_DIR!r} in the script, "
              f"resolved relative to the script's own folder, NOT your "
              f"terminal's current directory.)")
        print()
        print("  Files actually present in that folder:")
        try:
            entries = os.listdir(input_dir_resolved)
            if entries:
                for e in sorted(entries):
                    print(f"    {e}")
            else:
                print("    (folder is empty)")
        except FileNotFoundError:
            print(f"    (folder does not exist: {input_dir_resolved})")
        print()
        print("  Fix: edit INPUT_DIR near the top of this script to point at "
              "wherever your dtu_mhcii_batch_*.xls files actually live, "
              "relative to this script's folder.")
        sys.exit(1)

    print(f"Found {len(filepaths)} MHC-II batch files to parse...")

    all_dfs = []
    allele_panels_seen = {}  # alleles-string -> list of files, to catch inconsistent panels
    failures = []

    for fp in filepaths:
        name = os.path.basename(fp)
        print(f"  Processing {name}...")
        try:
            df = parse_dtu_mhcii_file(fp)
        except Exception as e:
            print(f"    ERROR: {e}")
            failures.append((name, str(e)))
            continue

        allele_key = df["_alleles"].iloc[0]
        allele_panels_seen.setdefault(allele_key, []).append(name)
        print(f"    OK -- {len(df)} peptides, alleles: {allele_key}")
        all_dfs.append(df)

    if not all_dfs:
        print("\nERROR: No files could be parsed. Fix the errors above before continuing.")
        sys.exit(1)

    # Flag (don't silently ignore) if files used different allele panels
    if len(allele_panels_seen) > 1:
        print("\nWARNING: not all files used the same allele panel:")
        for panel, files in allele_panels_seen.items():
            print(f"  {panel}\n    -> {files}")
        print("Concatenating anyway, but per-allele columns will NOT line up "
              "across files with different panels. Check '_alleles' and "
              "'_source_file' columns in the output before trusting cross-file comparisons.\n")

    combined = pd.concat(all_dfs, ignore_index=True, sort=False)

    if "NB" not in combined.columns:
        print("ERROR: 'NB' column missing after parsing -- cannot filter binders.")
        sys.exit(1)

    combined["NB"] = pd.to_numeric(combined["NB"], errors="coerce")
    strong_binders = combined[combined["NB"] >= MIN_ALLELES_HIT].copy()
    strong_binders = strong_binders.sort_values("NB", ascending=False)

    print(f"\nTotal peptides parsed: {len(combined)}")
    print(f"Strong binders (NB >= {MIN_ALLELES_HIT}): {len(strong_binders)}")

    if failures:
        print(f"\n{len(failures)} file(s) FAILED to parse:")
        for name, err in failures:
            print(f"  {name}: {err}")

    if strong_binders.empty:
        print("\nNo strong binders found among successfully parsed files. "
              "This is a real result, not a parsing failure -- don't reflexively "
              "assume the script is broken again. Check MIN_ALLELES_HIT and "
              "RANK_STRONG_THRESHOLD match what you actually want.")
    else:
        strong_binders.to_csv(OUTPUT_CSV, index=False)
        print(f"Strong binders written to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()