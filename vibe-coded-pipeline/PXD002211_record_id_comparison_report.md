# PXD002211: ID-free vs ionbot quantification, linked by record id

Date: 2026-10-29. Code: `compare_by_record_id.py` (this folder).

## Inputs

| Name | Producer | FlashLFQ | Rows |
|---|---|---|---|
| pipeline | `quantify_precursors.py`, server run of 2026-10-29 (`../PXD002211_results_20ppm/clanker-results/`), `--ppm 20 --chg --int` | 2.1.5 | 163,789 MS2 scans, 143,578 quantified |
| v4 | `mgf2FlashLFQ.py` (`QuantifiedPeaks_v4.tsv`), `--ppm 20 --chg --int` | 1.2.4 | 125,203 peaks, 163,789 scans |
| ionbot | `ionbot2FlashLFQ-comparison.py` (`QuantifiedPeaks.ionbot.comparison.tsv`), PSMs at 1% FDR, `--ppm 20 --chg --int` | 1.2.4 | 30,655 peaks, 47,419 PSMs |

All three files carry one record id per MS2 scan: `<sha1(run name)[:12]>_SCAN<scan>`.

## Record ids

- Every record id parses, and no scan is in more than one peak.
- In all three files the prefix equals `sha1(run name)` for 100% of rows.
- All 47,419 ionbot record ids occur, as exact strings, in the pipeline output and in v4. The
  pipeline and v4 contain the same 163,789 record ids.
- The comparison joins on (File Name, scan number), parsed from the record id. This gives the
  same result as a join on the full record id.

## Agreement

Intensity agreement is given as r² of log10 intensities: Pearson r², and r² against y = x.
"Same apex" means the same apex RT and the same apex m/z (within 1e-4).

| Comparison (per scan) | Both quantified | Same apex | r², all | r², same apex |
|---|---|---|---|---|
| pipeline vs v4 | 142,705 | 88.1% | 0.970 / 0.969 | 0.998 / 0.998 |
| ionbot vs pipeline | 40,827 | 69.1% | 0.623 / 0.542 | 0.998 / 0.998 |
| ionbot vs v4 | 40,854 | 68.6% | 0.629 / 0.556 | 0.999 / 0.999 |

- Pipeline vs v4: same numbers as the title-based comparison in
  `PXD002211_v3_comparison_report.md`, so the record ids link the scans correctly.
- Of the 43,498 PSMs that ionbot quantifies, the pipeline quantifies 40,827 (93.9%). v4 quantifies
  40,854. The pipeline and v4 perform almost the same against ionbot, so the FlashLFQ version
  (1.2.4 vs 2.1.5) is not the cause of the disagreement.
- When the apex is the same, the intensities are the same (r² 0.998).

## Cause of the disagreement: precursor isotope errors

FlashLFQ looks for the envelope at the mass that it is given. ionbot gives the peptide mass.
The ID-free runs give the precursor mass from the mzML selected-ion m/z. The offset between the
two, in isotopes (13C spacing), shows the problem:

| ionbot PSMs | 0 | -1 | -2 | +1 | other |
|---|---|---|---|---|---|
| all (47,419) | 68.0% | 12.4% | 5.2% | 4.2% | 10.3% |
| Thermo monoisotopic m/z present (29,798) | 84.4% | 2.5% | 3.1% | 3.7% | 6.2% |
| Thermo monoisotopic m/z absent (17,621) | 40.2% | 29.1% | 8.7% | 4.9% | 17.1% |

- The Thermo "Monoisotopic M/Z" trailer value is absent for 49% of the scans (40–65% per run).
  For those scans the mzML m/z is the isolation target, which is often not the monoisotopic peak.
  Offset -1 means that the precursor mass is one isotope too high.
- Offsets of 4 or more isotopes (7.5% of PSMs) are probably ionbot identifications of a
  co-isolated precursor. They were not examined.

The outcome for the 43,498 PSMs that ionbot quantifies:

| Offset | PSMs | Same apex | Other apex | Pipeline not quantified |
|---|---|---|---|---|
| 0 | 31,965 | 28,223 (88.3%) | 3,680 | 62 |
| not 0 | 11,533 | 1 | 8,923 | 2,609 |

- With offset 0, the pipeline and ionbot agree: r² 0.943 for all scans that both quantify.
- With an offset, the pipeline quantifies a different envelope or nothing. For these scans the
  r² is 0.065, and the pipeline intensity is a median 4.2× lower (log10 ratio +0.62).
- Offsets explain 2,609 of the 2,671 PSMs that the pipeline does not quantify, and 8,923 of the
  12,603 PSMs with a different apex. This answers the open question in the v3 report:
  the missing ionbot apexes come mainly from precursor isotope errors.
- 3,680 PSMs with offset 0 still have a different apex. Most have the same m/z and a different
  apex RT. The other likely contributor is the isotope model (peptide formula for ionbot,
  averagine for the pipeline). These cases were not examined further.

## Consequence

This is a limitation of the ID-free approach in its current form, not a FlashLFQ setting.
About 32% of the identified PSMs, and probably a similar share of all MS2 scans, are quantified at
the wrong mass. A monoisotopic precursor correction from the MS1 spectrum (find the first
isotope of the envelope before the isolation target) before the FlashLFQ step would address it.
This correction is not implemented.

## Reproduce

The three inputs are in `../PXD002211_results_20ppm/`. From this folder:

```bash
D=../PXD002211_results_20ppm
python compare_by_record_id.py --results $D/clanker-results \
    --v4 $D/QuantifiedPeaks_v4.tsv \
    --ionbot $D/QuantifiedPeaks.ionbot.comparison.tsv
```
