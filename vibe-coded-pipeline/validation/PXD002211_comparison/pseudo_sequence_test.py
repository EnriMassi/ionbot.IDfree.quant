"""Test whether the pseudo-sequence text changes FlashLFQ results.

Step 1 (make-input): take the FlashLFQ input that quantify_precursors.py wrote
for one run and replace each SYN_ id by the mgf2FlashLFQ.py pseudo-sequence
(scan number digits 0-9 -> A C D E F G H I L M). Masses, charges and RTs stay
identical, so only the sequence text differs.

Step 2: run FlashLFQ on that input with the same switches as the pipeline run,
for example:
    dotnet CMD.dll --idt v3style_input.tsv --rep <folder with the run's mzML> \
        --out out --ppm 20 --iso 5 --nis 2 --pipfdr 0 --rea --chg --int

Step 3 (compare): compare that QuantifiedPeaks.tsv with the pipeline's
QuantifiedPeaks.tsv for the same run and with QuantifiedPeaks_v3.tsv.
"""

import argparse
import os

import numpy as np
import pandas as pd

DIGIT_LETTERS = "ACDEFGHILM"  # same mapping as digit_to_aminoacid in mgf2FlashLFQ.py


def scan_to_sequence(scan: str) -> str:
    return "".join(DIGIT_LETTERS[int(d)] for d in scan)


def make_input(args: argparse.Namespace) -> None:
    run_input = os.path.join(args.results, "flashlfq", args.run, "flashlfq_input_run.tsv")
    records = pd.read_csv(run_input, sep="\t", dtype=str)
    inventory = pd.read_csv(os.path.join(args.results, "precursor_inventory.tsv"), sep="\t", dtype=str,
                            usecols=["synthetic_record_id", "scan_number"])
    merged = records.merge(inventory, left_on="Full Sequence", right_on="synthetic_record_id")
    assert len(merged) == len(records)
    merged["Full Sequence"] = merged["scan_number"].map(scan_to_sequence)
    merged["Base Sequence"] = merged["Full Sequence"]
    assert merged["Full Sequence"].is_unique
    merged[records.columns].to_csv(args.out_input, sep="\t", index=False)
    merged[["synthetic_record_id", "Full Sequence"]].to_csv(args.out_map, sep="\t", index=False)
    print(f"wrote {len(merged)} records to {args.out_input}")


def read_peaks(path: str) -> pd.DataFrame:
    peaks = pd.read_csv(path, sep="\t", low_memory=False)
    for column in ["Peak intensity", "Peak RT Apex", "Peak MZ", "PSMs Mapped", "Precursor Charge"]:
        peaks[column] = pd.to_numeric(peaks[column], errors="coerce")
    return peaks


def per_title(peaks: pd.DataFrame, ids: pd.Series) -> pd.DataFrame:
    out = peaks.assign(t=ids).explode("t")
    return out[["t", "Peak intensity", "Peak RT Apex", "Peak MZ"]]


def agreement(x: pd.DataFrame, y: pd.DataFrame, name: str) -> None:
    m = x.merge(y, on="t", suffixes=("_x", "_y"))
    q = (m["Peak intensity_x"] > 0) & (m["Peak intensity_y"] > 0)
    same = q & np.isclose(m["Peak RT Apex_x"], m["Peak RT Apex_y"], atol=1e-4) \
        & np.isclose(m["Peak MZ_x"], m["Peak MZ_y"], atol=1e-4)
    equal = same & np.isclose(m["Peak intensity_x"], m["Peak intensity_y"], rtol=1e-6)
    print(f"{name}: spectra {len(m)}, quantified in both {q.sum()}, "
          f"same apex {same.sum() / q.sum():.3f}, identical intensity {equal.sum() / q.sum():.3f}")


def compare(args: argparse.Namespace) -> None:
    seq_to_syn = pd.read_csv(args.map, sep="\t").set_index("Full Sequence")["synthetic_record_id"]
    syn_to_seq = pd.Series(seq_to_syn.index, index=seq_to_syn.values)
    title = pd.read_csv(os.path.join(args.results, "precursor_inventory.tsv"), sep="\t",
                        usecols=["synthetic_record_id", "spectrum_title"]).set_index("synthetic_record_id")["spectrum_title"]

    digit = read_peaks(args.digit_peaks)
    syn = read_peaks(os.path.join(args.results, "flashlfq", args.run, "output", "QuantifiedPeaks.tsv"))
    v3 = read_peaks(args.v3)
    v3 = v3[v3["File Name"] == args.run]

    split = lambda s: s.astype(str).str.split("|")
    digit_t = per_title(digit, split(digit["Full Sequence"]).map(lambda ids: [title[seq_to_syn[i]] for i in ids]))
    syn_t = per_title(syn, split(syn["Full Sequence"]).map(lambda ids: [title[i] for i in ids]))
    v3_t = per_title(v3, v3["Protein Group"].astype(str).str.split(";"))

    agreement(digit_t, syn_t, "digit letters (2.1.5) vs SYN ids (2.1.5)")
    agreement(digit_t, v3_t, "digit letters (2.1.5) vs v3 (1.2.4)    ")
    agreement(syn_t, v3_t, "SYN ids (2.1.5)       vs v3 (1.2.4)    ")

    # Dose response: single-spectrum peaks, by number of sulfur letters (C, M).
    d1 = digit[(digit["PSMs Mapped"] == 1) & (digit["Peak intensity"] > 0)]
    s1 = syn[(syn["PSMs Mapped"] == 1) & (syn["Peak intensity"] > 0)].copy()
    s1["seq"] = s1["Full Sequence"].map(syn_to_seq)
    m = d1.merge(s1, left_on="Full Sequence", right_on="seq", suffixes=("_d", "_s"))
    same = np.isclose(m["Peak RT Apex_d"], m["Peak RT Apex_s"], atol=1e-4) \
        & np.isclose(m["Peak MZ_d"], m["Peak MZ_s"], atol=1e-4)
    table = pd.DataFrame({
        "sulfur_letters": m["seq"].str.count("[CM]"),
        "log2_digit_over_syn": np.log2(m["Peak intensity_d"] / m["Peak intensity_s"]),
        "same_apex": same,
    })
    print(f"\nsingle-spectrum peaks in both: {len(m)}")
    summary = table.groupby("sulfur_letters").agg(
        n=("same_apex", "size"),
        same_apex=("same_apex", "mean"),
        median_log2_ratio_same_apex=("log2_digit_over_syn", lambda r: r[table.loc[r.index, "same_apex"]].median()),
    )
    print(summary.round(3))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    mk = sub.add_parser("make-input")
    mk.add_argument("--results", required=True, help="quantify_precursors.py output folder (--keep-intermediate not needed)")
    mk.add_argument("--run", required=True, help="run name, e.g. HUVEC_IL1b1b_sn")
    mk.add_argument("--out-input", default="v3style_input.tsv")
    mk.add_argument("--out-map", default="v3style_map.tsv")
    cp = sub.add_parser("compare")
    cp.add_argument("--results", required=True)
    cp.add_argument("--run", required=True)
    cp.add_argument("--digit-peaks", required=True, help="QuantifiedPeaks.tsv from the step 2 FlashLFQ run")
    cp.add_argument("--map", default="v3style_map.tsv")
    cp.add_argument("--v3", required=True, help="QuantifiedPeaks_v3.tsv")
    args = parser.parse_args()
    make_input(args) if args.command == "make-input" else compare(args)


if __name__ == "__main__":
    main()
