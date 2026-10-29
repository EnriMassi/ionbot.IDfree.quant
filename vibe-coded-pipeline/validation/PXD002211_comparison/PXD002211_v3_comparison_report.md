# PXD002211: old ID-free quantification (v3) vs quantify_precursors.py

Date: 2026-10-28. Data: the 6 HCD runs of PXD002211 (HUVEC_IL1b1a_sn … HUVEC_IL1b3b_sn).

## Runs compared

| Name | Producer | FlashLFQ | Switches | Pseudo-sequence |
|---|---|---|---|---|
| **v3** | `mgf2FlashLFQ.py` → `QuantifiedPeaks_v3.tsv` | 1.2.4 (Docker) | `--ppm 20 --chg --int` | scan digits as residues (`16080` → `CHALA`) |
| **pipeline** | `quantify_precursors.py` → `clanker-results/` | 2.1.5 (Docker, pinned digest) | `--ppm 20` (no `--chg`, no `--int`) | `SYN_` + SHA-1 hash (not valid residues) |
| **pipeline, v3 switches** | `quantify_precursors.py --ppm 20 --id-specific-charge --integrate` | 2.1.5 (local release build, same engine DLLs as the image) | `--ppm 20 --chg --int` | `SYN_` |

All runs use one FlashLFQ run per file, no MBR and no normalization. The spectra are joined on
spectrum title. v3 titles come from `Protein Group` (split on `;`). The pipeline has one row per spectrum.

## Agreement

Both cover the same 163,789 spectra. Every title is present in both, and every v3 title sits in exactly one peak.

| | pipeline (default) | pipeline, v3 switches |
|---|---|---|
| Same quantified / not-quantified call | 99.01% | 99.04% |
| Quantified in both | 142,662 | 142,705 |
| Quantified only in v3 / only in pipeline | 744 / 873 | 701 / 873 |
| Same apex centroid (RT and m/z), among quantified in both | 76.3% | 88.1% |
| log10(v3 / pipeline) intensity, same apex, median (IQR) | 0.91 (0.55–1.16) | 0.002 (0.000–0.015) |
| Pearson r, log10 intensity, same apex | 0.923 | 0.999 |
| Pearson r, log10 intensity, all quantified in both | 0.881 | 0.985 |
| **R², log10 intensity, all quantified in both** | **0.776** | **0.970** |
| R² against the line y = x (same scale, log10) | 0.296 | 0.969 |
| Feature charge ≠ precursor charge | 13.4% | 0.03% |
| Spectra on a shared feature (v3: 57,481) | 71,199 | 57,887 |
| ionbot (ID-based) peaks whose apex is found (v3: 73.4%) | 57.1% | 69.0% |

Quantification rates per file agree within 0.3 percentage points (83–91%).
The v3 rows with intensity 0 (20,383) all have no apex. They are the same case as the
pipeline's `no_ms1_feature_found` (20,254), not measured zeros.

## Causes of disagreement

Ordered by size of effect.

1. **`--int` (intensity definition).** v3 integrates the isotope envelope over the whole
   chromatographic peak. The pipeline default reports the apex envelope only. This explains the
   ~8× (10^0.91) intensity difference and most of the lower correlation. With `--integrate`, the
   median ratio is 1.00.

2. **`--chg` (charge state).** Without `--chg`, FlashLFQ may build the feature on another charge
   state of the same peptide (13.4% of spectra). It then chooses a different apex and merges
   spectra of different charges into one feature (71,199 vs 57,481 spectra on shared features).
   With `--id-specific-charge`, same-apex agreement rises from 76% to 88%, and ionbot overlap rises
   from 57% to 69%.

3. **Pseudo-sequence letters (isotope model).** The v3 pseudo-sequences are real amino acids.
   FlashLFQ therefore computes the isotope envelope from that pseudo-peptide's formula and fills the
   remaining mass with averagine. The `SYN_` ids are not valid sequences, so FlashLFQ uses pure
   averagine scaled to the precursor mass. Digits 1 and 9 become C and M, which add sulfur, so
   the v3 envelope depends on the scan number. Test on HUVEC_IL1b1b_sn: FlashLFQ 2.1.5, identical
   input except the sequence text.

   | Comparison | Same apex | Identical intensity |
   |---|---|---|
   | scan-digit letters vs `SYN_` (both 2.1.5) | 90.8% | 0.6% |
   | scan-digit letters (2.1.5) vs v3 (1.2.4) | 98.0% | 85.2% |
   | `SYN_` (2.1.5) vs v3 (1.2.4) | 89.8% | 0.5% |

   The effect grows with the number of C/M letters. Single-spectrum peaks:

   | C/M letters | 0 | 1 | 2 | 3 | 4 |
   |---|---|---|---|---|---|
   | Same apex as `SYN_` | 94% | 93% | 89% | 84% | 75% |
   | Median log2(letters / `SYN_`) intensity, same apex | 0.00 | +0.003 | +0.028 | +0.068 | +0.10 |

   This is the main remaining difference once the switches match. The bias is small (up to ~7%
   intensity), but it is arbitrary. Averagine is the neutral model for a precursor of unknown
   composition.

4. **FlashLFQ version (1.2.4 vs 2.1.5).** This effect is minor. With identical pseudo-sequences,
   98.0% of spectra have the same apex, and 85.2% have a bit-identical intensity.

When the apex differs, most differences are a different apex scan at the same m/z, or a
neighbouring isotope peak (±1.003/z). The remaining 1% of status disagreements (about 1,600
spectra) were not examined one by one. They are probably borderline features, where only one
run finds the 2 required isotope peaks.

## Conclusions

- With matching switches, the two implementations agree closely: 99% the same quantified call,
  88% the same apex, and r = 0.999 on the same apex. The remaining differences come mainly from
  the isotope model that the v3 pseudo-sequence letters introduce.
- To compare with the ionbot ID-based run (FlashLFQ `--chg --int`), run the pipeline with
  `--id-specific-charge --integrate`. These became the pipeline defaults on 2026-10-29.
- Keep the `SYN_` ids (or any non-residue id). Ids made of amino-acid letters change the
  quantification.
- Why 27–31% of ionbot peaks have no matching ID-free apex in either implementation: mainly
  precursor isotope errors in the ID-free input. See `PXD002211_record_id_comparison_report.md`.

## Reproduce

```bash
D=/path/to/PXD002211
# 1. Agreement tables (both result folders, optional ionbot overlap)
python compare_v3_vs_pipeline.py --v3 $D/QuantifiedPeaks_v3.tsv \
    --results $D/clanker-results --results $D/results_chg_int \
    --ionbot $D/QuantifiedPeaks.ionbot.tsv

# results_chg_int was made with:
python quantify_precursors.py --project $D --output $D/results_chg_int \
    --ppm 20 --id-specific-charge --integrate

# 2. Pseudo-sequence test on one run
python pseudo_sequence_test.py make-input --results $D/results_chg_int --run HUVEC_IL1b1b_sn
mkdir rep out && gunzip -c $D/MZML/HUVEC_IL1b1b_sn.mzML.gz > rep/HUVEC_IL1b1b_sn.mzML
docker run --rm -v "$PWD":/w -w /w \
    smithchemwisc/flashlfq:2.1.5@sha256:5340791eba72ce6c2f5d82e8d6d43a99b6498eb4de9244c213566101448adb7a \
    --idt v3style_input.tsv --rep rep --out out \
    --ppm 20 --iso 5 --nis 2 --pipfdr 0 --rea --chg --int
python pseudo_sequence_test.py compare --results $D/results_chg_int --run HUVEC_IL1b1b_sn \
    --digit-peaks out/QuantifiedPeaks.tsv --v3 $D/QuantifiedPeaks_v3.tsv
```

The `results_chg_int` run and the pseudo-sequence test in this report were made on Windows with
the FlashLFQ 2.1.5 release build (`--flashlfq-local-cmd`). Its engine DLLs are byte-identical to
the Docker image.
