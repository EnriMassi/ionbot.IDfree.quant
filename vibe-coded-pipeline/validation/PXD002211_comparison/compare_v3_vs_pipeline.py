"""Compare the old identification-free FlashLFQ output (QuantifiedPeaks_v3.tsv,
from mgf2FlashLFQ.py) with the output of quantify_precursors.py, per spectrum.

Usage:
    python compare_v3_vs_pipeline.py --v3 QuantifiedPeaks_v3.tsv \
        --results clanker-results [--results other-results ...] \
        [--ionbot QuantifiedPeaks.ionbot.tsv]

The join key is the spectrum title. In v3 the titles of a peak are in
'Protein Group', joined by ';'. The pipeline has one row per spectrum.
"""

import argparse
import os

import numpy as np
import pandas as pd

NUMERIC_PEAK_COLUMNS = ["Peak intensity", "Peak RT Apex", "Peak MZ", "Peak Charge", "PSMs Mapped"]


def read_peaks(path: str) -> pd.DataFrame:
    peaks = pd.read_csv(path, sep="\t", low_memory=False)
    for column in NUMERIC_PEAK_COLUMNS:
        peaks[column] = pd.to_numeric(peaks[column], errors="coerce")
    return peaks


def v3_per_spectrum(v3: pd.DataFrame) -> pd.DataFrame:
    """One row per spectrum title. Pairs titles with the peak, never with the
    position of a pseudo-sequence in 'Full Sequence'."""
    titles = v3["Protein Group"].astype(str).str.split(";")
    n_seq = v3["Full Sequence"].astype(str).str.split("|").str.len()
    assert (titles.str.len() == v3["PSMs Mapped"]).all(), "titles != PSMs Mapped"
    assert (titles.str.len() == n_seq).all(), "titles != pseudo-sequences"
    spectra = v3.assign(spectrum_title=titles).explode("spectrum_title")
    assert spectra["spectrum_title"].is_unique, "a title appears in more than one v3 peak"
    return spectra[["spectrum_title", "Peak intensity", "Peak RT Apex", "Peak MZ", "Peak Charge", "PSMs Mapped"]]


def apex_key(run: pd.Series, rt: pd.Series, mz: pd.Series) -> pd.Series:
    """Apex centroid identity: same run, same apex scan RT, same apex m/z."""
    return run + "|" + rt.round(4).astype(str) + "|" + mz.round(4).astype(str)


def log10_pearson(a: pd.Series, b: pd.Series) -> float:
    return float(np.corrcoef(np.log10(a), np.log10(b))[0, 1])


def compare(v3: pd.DataFrame, v3_spectra: pd.DataFrame, results_dir: str, ionbot: pd.DataFrame | None) -> None:
    new = pd.read_csv(os.path.join(results_dir, "ms1_precursor_quantification.tsv"), sep="\t", low_memory=False)
    print(f"\n######## {results_dir}")
    print("pipeline status counts:", new["quantification_status"].value_counts().to_dict())

    joined = new.merge(v3_spectra, on="spectrum_title", how="outer", indicator=True)
    print("title join:", joined["_merge"].value_counts().to_dict())
    joined = joined[joined["_merge"] == "both"].copy()
    joined["old_q"] = joined["Peak RT Apex"].notna() & (joined["Peak intensity"] > 0)
    joined["new_q"] = joined["quantification_status"] == "quantified"

    print("\nquantified (rows: pipeline, columns: v3):")
    print(pd.crosstab(joined["new_q"], joined["old_q"]))
    agree = (joined["old_q"] == joined["new_q"]).mean()
    print(f"status agreement: {agree:.4f}")
    print("\nquantified fraction per file:")
    print(joined.groupby("source_file")[["old_q", "new_q"]].mean().round(3))

    both = joined[joined["old_q"] & joined["new_q"]].copy()
    same_apex = (
        np.isclose(both["Peak RT Apex"], both["feature_rt_apex"], atol=1e-4)
        & np.isclose(both["Peak MZ"], both["feature_mz"], atol=1e-4)
    )
    print(f"\nquantified in both: {len(both)}")
    print(f"same apex centroid: {same_apex.sum()} ({same_apex.mean():.3f})")
    shift = ((both["Peak MZ"] - both["feature_mz"]) * both["precursor_charge"])[~same_apex].round(2)
    print("apex differs, (v3 m/z - pipeline m/z) * z, most common:", shift.value_counts().head(6).to_dict())

    ratio = np.log10(both.loc[same_apex, "Peak intensity"] / both.loc[same_apex, "ms1_intensity"])
    print("same apex, log10(v3 / pipeline) intensity quantiles:",
          ratio.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).round(3).to_dict())
    print(f"pearson r log10 intensity, same apex: {log10_pearson(both.loc[same_apex, 'Peak intensity'], both.loc[same_apex, 'ms1_intensity']):.4f}")
    print(f"pearson r log10 intensity, all quantified in both: {log10_pearson(both['Peak intensity'], both['ms1_intensity']):.4f}")

    charge_diff = (both["feature_charge"] != both["precursor_charge"]).sum()
    print(f"pipeline feature charge != precursor charge: {charge_diff} ({charge_diff / len(both):.3f})")

    shared_new = new[new["feature_n_spectra"] > 1]
    print(f"\nspectra on a shared feature: pipeline {len(shared_new)} ({shared_new['feature_id'].nunique()} features, "
          f"max {int(new['feature_n_spectra'].max())}); v3 {(v3_spectra['PSMs Mapped'] > 1).sum()} "
          f"({(v3['PSMs Mapped'] > 1).sum()} peaks, max {int(v3['PSMs Mapped'].max())})")

    if ionbot is None:
        return
    ion = ionbot[ionbot["Peak RT Apex"].notna() & (ionbot["Peak intensity"] > 0)]
    ion_keys = apex_key(ion["File Name"], ion["Peak RT Apex"], ion["Peak MZ"])
    feats = new[new["quantification_status"] == "quantified"].drop_duplicates(["source_file", "feature_id"])
    feat_run = feats["source_file"].str.replace(r"\.mzML(\.gz)?$", "", regex=True)
    new_keys = set(apex_key(feat_run, feats["feature_rt_apex"], feats["feature_mz"]))
    v3q = v3[v3["Peak RT Apex"].notna() & (v3["Peak intensity"] > 0)]
    v3_keys = set(apex_key(v3q["File Name"], v3q["Peak RT Apex"], v3q["Peak MZ"]))
    print(f"\nionbot quantified peaks: {len(ion)}; apex found in pipeline: {ion_keys.isin(new_keys).mean():.3f}; "
          f"in v3: {ion_keys.isin(v3_keys).mean():.3f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--v3", required=True, help="QuantifiedPeaks_v3.tsv")
    parser.add_argument("--results", required=True, action="append", help="quantify_precursors.py output folder")
    parser.add_argument("--ionbot", help="identification-based QuantifiedPeaks.tsv (optional)")
    args = parser.parse_args()

    v3 = read_peaks(args.v3)
    print(f"v3 peaks: {len(v3)}; with apex: {v3['Peak RT Apex'].notna().sum()}; "
          f"intensity 0: {(v3['Peak intensity'] == 0).sum()}; "
          f"intensity 0 with apex: {((v3['Peak intensity'] == 0) & v3['Peak RT Apex'].notna()).sum()}")
    v3_spectra = v3_per_spectrum(v3)
    print(f"v3 spectra: {len(v3_spectra)}")
    ionbot = read_peaks(args.ionbot) if args.ionbot else None
    for results_dir in args.results:
        compare(v3, v3_spectra, results_dir, ionbot)


if __name__ == "__main__":
    main()
