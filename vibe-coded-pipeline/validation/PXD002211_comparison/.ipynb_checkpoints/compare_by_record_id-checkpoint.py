"""Compare three FlashLFQ quantifications of PXD002211 per MS2 scan, linked by
record id (<sha1(file name)[:12]>_SCAN<scan>):

- quantify_precursors.py output (ms1_precursor_quantification.tsv, one row per scan);
- mgf2FlashLFQ.py output (QuantifiedPeaks_v4.tsv, record ids '|'-joined in Full Sequence);
- ionbot2FlashLFQ-comparison.py output (QuantifiedPeaks.ionbot.comparison.tsv,
  one record id per ionbot PSM, '|'-joined in Full Sequence).

Usage:
    python compare_by_record_id.py --results clanker-results \
        --v4 QuantifiedPeaks_v4.tsv --ionbot QuantifiedPeaks.ionbot.comparison.tsv \
        [--ionbot-peptide QuantifiedPeaks.ionbot.tsv]

The join key is (File Name, scan number), parsed from the record id. The hash
prefix is only checked, not used, so files with a wrong prefix still join.
Earlier versions of ionbot2FlashLFQ-comparison.py hashed '<run>.mgf.gzip'.
"""

import argparse
import hashlib
import os

import numpy as np
import pandas as pd

C13 = 1.0033548378
PROTON = 1.007276466879
NUMERIC_PEAK_COLUMNS = ["Peak intensity", "Peak RT Apex", "Peak MZ", "Peak Charge",
                        "Peptide Monoisotopic Mass", "Precursor Charge"]


def file_hash(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


def peaks_per_scan(path: str, tag: str) -> pd.DataFrame:
    """One row per record id, with the values of the peak that contains it."""
    peaks = pd.read_csv(path, sep="\t", low_memory=False)
    for column in NUMERIC_PEAK_COLUMNS:
        peaks[column] = pd.to_numeric(peaks[column], errors="coerce")
    peaks["n_ids"] = peaks["Full Sequence"].astype(str).str.count(r"\|") + 1
    scans = peaks.assign(rid=peaks["Full Sequence"].astype(str).str.split("|")).explode("rid")
    scans["prefix"] = scans["rid"].str.extract(r"^([0-9a-f]{12})_SCAN")[0]
    scans["scan"] = pd.to_numeric(scans["rid"].str.extract(r"_SCAN(\d+)$")[0])
    scans["run"] = scans["File Name"]
    scans["q"] = scans["Peak RT Apex"].notna() & (scans["Peak intensity"] > 0)
    assert scans["scan"].notna().all(), f"{tag}: unparsed record ids"
    assert not scans.duplicated(["run", "scan"]).any(), f"{tag}: a scan is in more than one peak"
    print(f"[{tag}] peak rows {len(peaks)}, with apex {peaks['Peak RT Apex'].notna().sum()}, scans {len(scans)}")
    print(f"[{tag}] prefix == sha1(run): {(scans['prefix'] == scans['run'].map(file_hash)).mean():.3f}, "
          f"== sha1(run + '.mgf.gzip'): {(scans['prefix'] == (scans['run'] + '.mgf.gzip').map(file_hash)).mean():.3f}")
    return scans


def read_pipeline(results: str) -> pd.DataFrame:
    rows = pd.read_csv(os.path.join(results, "ms1_precursor_quantification.tsv"), sep="\t", low_memory=False)
    rows["run"] = rows["source_file"].str.replace(r"\.mzML(\.gz)?$", "", regex=True)
    rows["scan"] = pd.to_numeric(rows["scan_number"])
    rows["q"] = rows["quantification_status"] == "quantified"
    print(f"[pipeline] rows {len(rows)}, {rows['quantification_status'].value_counts().to_dict()}")
    print(f"[pipeline] prefix == sha1(run): {(rows['synthetic_record_id'].str[:12] == rows['run'].map(file_hash)).mean():.3f}")
    return rows


def r2(a: pd.Series, b: pd.Series) -> str:
    """Pearson r^2 and r^2 against y = x of log10 intensities, and median log10(a/b)."""
    a, b = np.log10(a.to_numpy()), np.log10(b.to_numpy())
    pearson = np.corrcoef(a, b)[0, 1] ** 2
    identity = 1 - ((a - b) ** 2).sum() / ((a - a.mean()) ** 2).sum()
    return f"n {len(a)}, r2 {pearson:.3f}, r2 vs y=x {identity:.3f}, median log10 ratio {np.median(a - b):+.3f}"


def same_apex(df: pd.DataFrame, a: tuple[str, str], b: tuple[str, str]) -> pd.Series:
    return (np.isclose(df[a[0]], df[b[0]], atol=1e-4) & np.isclose(df[a[1]], df[b[1]], atol=1e-4))


def compare(left: pd.DataFrame, right: pd.DataFrame, names: tuple[str, str],
            intensity: tuple[str, str], apex: tuple[tuple[str, str], tuple[str, str]]) -> pd.DataFrame:
    a, b = names
    merged = left.merge(right, on=["run", "scan"], how="outer", suffixes=(f"_{a}", f"_{b}"), indicator=True)
    print(f"\n== {a} vs {b}, per scan ==")
    print("join", merged["_merge"].value_counts().to_dict())
    merged = merged[merged["_merge"] == "both"].copy()
    print(pd.crosstab(merged[f"q_{a}"].rename(f"{a} quantified"), merged[f"q_{b}"].rename(f"{b} quantified")))
    both = merged[merged[f"q_{a}"] & merged[f"q_{b}"]]
    same = same_apex(both, *apex)
    print(f"both quantified {len(both)}, same apex (RT and m/z) {same.sum()} ({same.mean():.0%})")
    print("  all both quantified:", r2(both[intensity[0]], both[intensity[1]]))
    print("  same apex          :", r2(both.loc[same, intensity[0]], both.loc[same, intensity[1]]))
    return merged


def isotope_offsets(merged: pd.DataFrame) -> None:
    """Offset, in isotopes, between the ionbot peptide mass and the precursor mass
    that the ID-free run sends to FlashLFQ (from the mzML selected-ion m/z)."""
    m = merged.copy()
    m["iso"] = ((m["Peptide Monoisotopic Mass"] - m["neutral_mass"]) / C13).round().clip(-4, 4)
    has_mono = m["thermo_monoisotopic_mz"].notna()
    print("\n== ionbot peptide mass - ID-free precursor mass, in isotopes (|offset| >= 4 binned as +-4) ==")
    print("all PSMs              ", m["iso"].value_counts(normalize=True).sort_index().round(3).to_dict())
    print(f"Thermo mono m/z present ({has_mono.sum()})", m.loc[has_mono, "iso"].value_counts(normalize=True).sort_index().round(3).to_dict())
    print(f"Thermo mono m/z absent  ({(~has_mono).sum()})", m.loc[~has_mono, "iso"].value_counts(normalize=True).sort_index().round(3).to_dict())

    q = m[m["q_ionbot"]]
    group = np.where(~q["q_pipeline"], "pipeline not quantified",
                     np.where(same_apex(q, ("Peak RT Apex", "Peak MZ"), ("feature_rt_apex", "feature_mz")),
                              "same apex", "other apex"))
    print("\nionbot-quantified PSMs by isotope offset and outcome in the ID-free run:")
    print(pd.crosstab(q["iso"].rename("offset"), group, margins=True))
    for label, sel in (("offset 0", q["iso"] == 0), ("offset != 0", q["iso"] != 0)):
        sel = sel & q["q_pipeline"]
        print(f"  {label}, both quantified:", r2(q.loc[sel, "Peak intensity"], q.loc[sel, "ms1_intensity"]))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", required=True, help="quantify_precursors.py output folder")
    p.add_argument("--v4", required=True, help="QuantifiedPeaks_v4.tsv (mgf2FlashLFQ.py)")
    p.add_argument("--ionbot", required=True, help="QuantifiedPeaks.ionbot.comparison.tsv")
    p.add_argument("--ionbot-peptide", help="QuantifiedPeaks.ionbot.tsv (peptide-keyed run), to check the peaks are the same")
    args = p.parse_args()

    pipeline = read_pipeline(args.results)
    v4 = peaks_per_scan(args.v4, "v4")
    ionbot = peaks_per_scan(args.ionbot, "ionbot")

    if args.ionbot_peptide:
        old = pd.read_csv(args.ionbot_peptide, sep="\t", low_memory=False)
        key = lambda d: (d["File Name"] + "|" + pd.to_numeric(d["Peak RT Apex"]).round(4).astype(str)
                         + "|" + pd.to_numeric(d["Peak MZ"]).round(4).astype(str))
        k_old = set(key(old[pd.to_numeric(old["Peak RT Apex"], errors="coerce").notna()]))
        k_new = set(key(ionbot[ionbot["Peak RT Apex"].notna()].drop_duplicates("Full Sequence")))
        print(f"\n[ionbot peptide-keyed vs PSM-keyed] apex peaks {len(k_old)} vs {len(k_new)}, shared {len(k_old & k_new)}")

    peak_cols = ["run", "scan", "q", "Peak intensity", "Peak RT Apex", "Peak MZ"]
    pipe_cols = ["run", "scan", "q", "ms1_intensity", "feature_rt_apex", "feature_mz",
                 "neutral_mass", "thermo_monoisotopic_mz"]
    compare(pipeline[pipe_cols], v4[peak_cols], ("pipeline", "v4"), ("ms1_intensity", "Peak intensity"),
            (("feature_rt_apex", "feature_mz"), ("Peak RT Apex", "Peak MZ")))
    ion_cols = peak_cols + ["Peptide Monoisotopic Mass"]
    merged = compare(ionbot[ion_cols], pipeline[pipe_cols], ("ionbot", "pipeline"), ("Peak intensity", "ms1_intensity"),
                     (("Peak RT Apex", "Peak MZ"), ("feature_rt_apex", "feature_mz")))
    compare(ionbot[ion_cols], v4[peak_cols], ("ionbot", "v4"), ("Peak intensity_ionbot", "Peak intensity_v4"),
            (("Peak RT Apex_ionbot", "Peak MZ_ionbot"), ("Peak RT Apex_v4", "Peak MZ_v4")))
    isotope_offsets(merged)


if __name__ == "__main__":
    main()
