# Validation report: PXD009327

Date: 2026-10-29. Pipeline version 1.0.0, current defaults: `--id-specific-charge --integrate`
(FlashLFQ `--chg --int`), `--ppm 10`. Record IDs: `<SHA-1 of run name, 12 hex>_SCAN<scan>`.
Input: the three files in `PXD009327/MZML/`, run with `--project PXD009327`.
The full summary is in `PXD009327_quantification_summary.json`, and example rows are in
`example_ms1_precursor_quantification.tsv`.

This report replaces the version of 2026-10-28. That version used FlashLFQ's own defaults
(all charge states, apex intensity) and `SYN_<hash>` record IDs.

## What was and was not executed

| Component | Status |
|---|---|
| mzML extraction, inventory, FlashLFQ input generation, mapping, outputs | Executed on all 3 files |
| FlashLFQ quantification | Executed with the **official FlashLFQ 2.1.5 release build** (`dotnet CMD.dll`, .NET 10.0.14) through the pipeline's development runner (`--flashlfq-local-cmd`), because no Docker daemon runs on the development machine (Windows 11) |
| FlashLFQ inside Docker | Not executed on PXD009327. Docker execution was confirmed on a Linux server with PXD002211 (6 runs, 163,789 spectra, all checks passed, 2026-10-28). Run the server checklist below to confirm PXD009327 in Docker. |

The development run and the Docker image run the same FlashLFQ code. This was checked
without a Docker daemon, by streaming the image layers from Docker Hub:

- The image is `smithchemwisc/flashlfq:2.1.5@sha256:5340791e…adb7a` (linux/amd64, Alpine,
  .NET 10.0.12), with entrypoint `["dotnet", "/flashlfq/CMD.dll"]`.
- 118 of the 127 files in `/flashlfq` are byte-identical to the GitHub release 2.1.5 zip.
  These include `FlashLFQ.dll`, `MassSpectrometry.dll`, `Readers.dll` and `Chemistry.dll`,
  which hold the whole quantification engine (mzLib 1.0.592).
- `CMD.dll` and `Util.dll` (the CLI wrapper) are the same size and were built from the same git
  commit (`754dbd9859ae…`). Only the version stamp differs: `1.0.0.0+754dbd98` in the image
  and `2.1.5+754dbd98` in the release. The other differing files are the matching
  `.deps.json`, `.pdb`, `.config` and `README.txt` files.
- The local `FlashLFQ-master.zip` is the same commit, `754dbd98` (tag 2.1.5).

All 18 pipeline checks passed (exit code 0).

## Spectrum inventory (checks 1–4, 13)

| File | Spectra | MS1 | MS2 | HCD MS2 | PRIDE metadata MS2 | Record ID prefix |
|---|---|---|---|---|---|---|
| VEL180306-01_FWeiland180306H01_FW_180307072156 | 12311 | 2073 | 10238 | 10238 | 10238 | `7a46706c58fa_` |
| VEL180306-03_Fweiland180306H02_FW | 11905 | 2351 | 9554 | 9554 | n/a | `458053a626eb_` |
| VEL180306-05_Fweiland180306H03_FW | 11946 | 2330 | 9616 | 9616 | n/a | `d7190da27fc1_` |
| **Total** | | | | **29408** | | |

- All MS2 spectra are HCD (`MS:1000422`, filter `@hcd35.00`). Each has exactly one precursor and
  one selected ion, with charge and intensity, so all 29408 have precursor status `ok`.
- RT is given in minutes (`UO:0000031`), so no conversion was needed.
- The selected ion m/z equals Thermo's "Monoisotopic M/Z" trailer for every spectrum (0 differ by
  more than 0.01). The mzML files were written by ThermoRawFileParser 1.4.5.
- The mzML files contain **no spectrum title** (`MS:1000796`). Titles were therefore rebuilt as
  `"<mzML run id>: <nativeID>"`. All 29408 rebuilt titles match the `TITLE=` lines of the PRIDE
  MGF files in `PXD009327/MGF` exactly, in both directions. For every spectrum the MGF also
  agrees on PEPMASS (< 1e-4), CHARGE and RTINSECONDS (< 0.001 min).
- Internal IDs (`file::nativeID`), titles and record IDs are all unique. Every file has unique
  scan numbers, so no file needed the spectrum-index fallback. All 29408 record IDs were
  regenerated independently from `source_file` and `scan_number` with
  `sha1(run name)[:12] + "_SCAN" + scan`, and all matched.

## FlashLFQ input and reader (checks 5, 6)

- 29408 records were written (one per eligible spectrum), each with `Base Sequence` =
  `NO_SEQUENCE` (averagine isotope model) and the record ID as `Full Sequence`. They passed the
  pipeline's re-implementation of the Generic-reader constraints.
- FlashLFQ reported `Done reading PSMs; found N` with N = 10238, 9554 and 9616, which equals the
  records sent. There were no reader rejection messages.
- The effective `FlashLfqSettings.toml` was verified for each run: Normalize=false,
  Integrate=true, IdSpecificChargeState=true, PpmTolerance=10, IsotopePpmTolerance=5,
  NumIsotopesRequired=2, BayesianProteinQuant=false and MbrQValueThreshold=0.
- No MBR peaks were produced. Every QuantifiedPeaks row has detection type `MSMS`.
- **Control experiment** (602 records, 301 each from runs 01 and 03, with the current record
  IDs and `--chg --int`):
  - With stock MBR settings, one joint FlashLFQ run produced **67 MBR peaks**. These are
    synthetic records transferred into the other run.
  - With `--pipfdr 0` it produced none. Its 568 MSMS peak rows are identical to those of the
    stock run.
  - The joint `--pipfdr 0` run and two per-file runs gave identical MSMS peak tables
    (568 of 568 rows).

## Quantification and mapping (checks 7–12)

| | Count |
|---|---|
| Output rows (= HCD spectra) | 29408 (10238 / 9554 / 9616 per file, equal to the inventory) |
| `quantified` | 25049 (85.2%; 85.4 / 85.7 / 84.5% per file) |
| `no_ms1_feature_found` | 4359 |
| Failed or missing from FlashLFQ output | 0 |
| FlashLFQ peak rows | 25816, of which 21457 have an apex (= features) |
| Features with more than 1 spectrum | 2395 features covering 5987 spectra (maximum 15 spectra) |
| Mapping problems (unknown, duplicate or cross-file IDs, PSMs Mapped mismatches) | 0 |
| Single-record peaks echoing input mass, charge and RT | All |

Quantification rate by precursor charge: z2 87.1%, z3 84.8%, z4 85.4%, z5 84.9%, z6 77.2%,
z7 72.5%, z8 73.4%.

**Independent check against the raw MS1 data.** 300 random features from run 03 were checked
against the raw MS1 scans:
- **Apex peak.** The FlashLFQ apex m/z was present in the raw MS1 scan at the apex RT within
  0.5 ppm for 300 of 300 features.
- **Intensity.** For each MS1 scan between the feature's RT start and end, the observed isotope
  peaks (monoisotopic to 3 isotopes above the apex isotope, 5 ppm) were summed and divided by z.
  The ratio of that integrated sum to the FlashLFQ intensity had median 0.97 (IQR 0.71–1.05).

The intensity check is an approximation. FlashLFQ uses only the isotopes its averagine model
expects, imputes isotopes that are missing, and skips scans whose envelope fails its
correlation test. The ratio therefore spreads more than the apex-only check of the previous
version (median 1.02, IQR 0.97–1.06). The median confirms that `--int` reports the envelope
intensity (isotope sum / z) summed over the peak's MS1 scans (see README).

**Example of a shared feature** (run 01). Three 4+ spectra of the same precursor, taken over
3.5 minutes, were mapped to one FlashLFQ feature:

| scan | record ID | precursor m/z | z | MS2 RT | feature | feature m/z | RT start–apex–end | intensity |
|---|---|---|---|---|---|---|---|---|
| 6498 | `7a46706c58fa_SCAN6498` | 695.8583 | 4 | 64.592 | …#F003779 | 696.1092 | 62.67–63.10–68.30 | 3.58e6 |
| 6618 | `7a46706c58fa_SCAN6618` | 695.8578 | 4 | 65.695 | …#F003779 | 696.1092 | 62.67–63.10–68.30 | 3.58e6 |
| 6890 | `7a46706c58fa_SCAN6890` | 695.8544 | 4 | 68.131 | …#F003779 | 696.1092 | 62.67–63.10–68.30 | 3.58e6 |

The feature m/z is the apex isotope peak, here the M+1 isotope (+1.0034/4), not the
monoisotopic peak. FlashLFQ traces the most abundant isotope.

## Attribution diagnostics (reported, not filtered)

- **Charge state.** With `--chg`, only the precursor's own charge is searched. 18 of 25049
  quantified spectra (0.07%) have a feature charge that differs from their precursor charge. All
  18 are on shared features: two spectra at the same m/z received different charge
  assignments from the instrument, so their records reached the same apex peak and FlashLFQ
  merged them. The feature takes the charge of one of the merged records. (FlashLFQ's own default,
  without `--chg`, gave 12.8% here.)
- **RT bounds.** For 6510 quantified spectra (26.0%), the MS2 RT lies outside FlashLFQ's peak
  bounds. The median |MS2 RT − apex| is 0.14 min (IQR 0.06–0.50 min). The median feature width
  (end − start) is 0.65 min (IQR 0.19–2.29). A long tail comes from persistent ions; the largest
  feature (15 spectra, 6+, m/z 465.596) spans 82.4–102.7 min. The
  `ms2_rt_within_feature_bounds` and `ms2_rt_minus_feature_apex` columns let you filter these.
- **Apex mass error.** Median |error| is 0.92 ppm, and the maximum is 10 ppm (the tolerance).
- **Intensity range.** Feature intensities span log10 2.3–8.9 (median 4.6).

## Reproducibility

Two runs, with 16 and 4 threads, produced byte-identical `ms1_precursor_quantification.tsv`,
`feature_mapping.tsv`, `precursor_inventory.tsv` and `flashlfq_input.tsv`. Run time on a
16-core workstation was 47 s in total (extraction 4.7 s, FlashLFQ about 11 s per run).

## Automated tests

`python -m pytest tests` gives 22 passed. The tests cover:

- mass conversion
- record IDs (file hash + scan, index fallback) and the `NO_SEQUENCE` base sequence
- `--chg` and `--int` on by default, and the `--no-` switches
- HCD selection, CID and EThcD exclusion, ambiguous precursors, and RT unit conversion
- title preservation and reconstruction
- duplicate IDs
- Generic-reader constraints
- the bool-switch-safe CLI
- Docker mounts and quoting (spaces and commas)
- merged-peak mapping and integrity errors
- status distinctions
- end-to-end orchestration with a fake FlashLFQ, including silent failures and timeouts
- `--dry-run`
- `--project` path derivation
- PXD009327 HCD counts

## Server-side checklist (Docker)

```bash
docker pull smithchemwisc/flashlfq:2.1.5@sha256:5340791eba72ce6c2f5d82e8d6d43a99b6498eb4de9244c213566101448adb7a
python3 -m pytest tests                       # expect 22 passed
python3 quantify_precursors.py --project /data/PXD009327 --output /data/PXD009327/results_docker \
    --threads 16 --ppm 10                     # MZML/ input and MGF/ title check come from --project
echo $?                                       # expect 0 (2 = a validation check failed, 1 = fatal error)
python3 -c "import json;s=json.load(open('/data/PXD009327/results_docker/quantification_summary.json'));print(s['all_checks_passed'], s['counts']['quantification_status'])"
# expect: True {'quantified': 25049, 'no_ms1_feature_found': 4359}
# Optional: compare with a development-run table
cmp results_dev/ms1_precursor_quantification.tsv /data/PXD009327/results_docker/ms1_precursor_quantification.tsv
ls -l /data/PXD009327/results_docker          # files must be owned by the invoking user
```

The development run was on Windows and the server is Linux. The Linux .NET runtime should
produce the same counts. If the counts differ, compare the per-run `QuantifiedPeaks.tsv` files;
any difference would come from the platform, not from the pipeline.
