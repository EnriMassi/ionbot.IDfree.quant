# Identification-free MS1 precursor quantification with FlashLFQ

`quantify_precursors.py` measures the MS1 chromatographic feature behind **every HCD MS2
spectrum** in a set of mzML runs. It does not need peptide identifications. FlashLFQ, running
in Docker, does all peak detection and intensity extraction. Each spectrum's result is written
against its original spectrum title, so you can join it afterwards with PSMs from any search
engine.

```
mzML(.gz) ──► precursor_inventory.tsv ──► flashlfq_input.tsv (synthetic records)
          ──► FlashLFQ 2.1.5 in Docker, one container per run ──► QuantifiedPeaks.tsv
          ──► identifier-based mapping ──► ms1_precursor_quantification.tsv (1 row per HCD spectrum)
```

## Quick start (Linux server)

Requirements:

- Python ≥ 3.11 (standard library only)
- Docker (the user must be allowed to run `docker`)
- About 2× the uncompressed mzML size of free space in the temporary directory, if the inputs
  are `.mzML.gz`

```bash
docker pull smithchemwisc/flashlfq:2.1.5@sha256:5340791eba72ce6c2f5d82e8d6d43a99b6498eb4de9244c213566101448adb7a

# Project folder with the standard layout: <project>/MZML, <project>/MGF
python3 quantify_precursors.py --project /data/PXD009327 --threads 16 --ppm 10
#   -> reads /data/PXD009327/MZML, writes /data/PXD009327/IDfree-quant-results,
#      checks titles against /data/PXD009327/MGF (if it contains .mgf/.mgf.gz files)

# Equivalent with explicit paths (any explicit path overrides the derived one)
python3 quantify_precursors.py --input /data/PXD009327/MZML --output /data/PXD009327/IDfree-quant-results \
    --reference-mgf-dir /data/PXD009327/MGF --threads 16 --ppm 10

python3 quantify_precursors.py --project /data/PXD009327 --dry-run   # print the plan and docker commands
python3 -m pytest tests                                              # 22 tests, no Docker needed
```

Exit codes:

- `0`: all checks passed
- `1`: fatal error (bad input, Docker unavailable, …)
- `2`: outputs were written, but at least one validation check failed or a run failed. See
  `quantification_summary.json` → `checks`.

### Options

| Option | Default | Meaning |
|---|---|---|
| `--project` | none | Project folder. Derives `--input <project>/MZML`, `--output <project>/IDfree-quant-results` and `--reference-mgf-dir <project>/MGF`. The MGF check is used only when that folder holds `.mgf`/`.mgf.gz` files. Folder names are matched exactly first, then case-insensitively. |
| `--input` | `<project>/MZML` | Directory with `.mzML` / `.mzML.gz` files. Only read, never modified. Required without `--project`. |
| `--output` | `<project>/IDfree-quant-results` | Results directory. Created if missing. Required without `--project`. |
| `--threads` | 4 | FlashLFQ `--thr`, Docker `--cpus`, and the number of parallel extraction processes |
| `--ppm` | 10 | FlashLFQ peak-finding mass tolerance (`--ppm`) |
| `--isotope-ppm` | 5 | FlashLFQ isotope tolerance (`--iso`) |
| `--num-isotopes` | 2 | Isotope peaks required (`--nis`, minimum 2) |
| `--id-specific-charge` / `--no-id-specific-charge` | on | Search only the precursor's own charge state (`--chg`). `--no-id-specific-charge` searches all charges, which is FlashLFQ's own default. |
| `--integrate` / `--no-integrate` | on | Sum the envelopes over the chromatographic peak (`--int`). `--no-integrate` reports the apex envelope only, which is FlashLFQ's own default. |
| `--image` | pinned 2.1.5 digest | FlashLFQ image |
| `--docker-memory` | none | Docker `--memory` limit, for example `32g` |
| `--docker-arg` | none | Extra `docker run` argument, repeatable. For example, on SELinux hosts use `--docker-arg=--security-opt=label=disable`. |
| `--timeout` | none | Seconds per FlashLFQ run. When it expires, the container is killed and the run is marked failed. |
| `--tmp-dir` | system temp | Where `.mzML.gz` files are decompressed for FlashLFQ |
| `--keep-intermediate` | off | Keep the decompressed staging files |
| `--title-template` | `{run_id}: {native_id}` | Used only when the mzML has no `MS:1000796` spectrum title (see below) |
| `--reference-mgf-dir` | `<project>/MGF` | Check the titles, m/z, charge and RT against the `TITLE=` lines of existing MGF files |
| `--no-reference-mgf` | off | Skip the MGF check, even when `<project>/MGF` exists |
| `--extract-only` | off | Stop after writing the inventory and the FlashLFQ input |
| `--dry-run` | off | Print the plan and the docker commands, then exit without reading or writing anything |
| `--log-level` | info | Console log level. The log files always record debug detail. |

## Outputs (`--output`)

| File | Content |
|---|---|
| `ms1_precursor_quantification.tsv` | **Main result.** One row per HCD MS2 spectrum, in inventory order, including spectra that failed. |
| `precursor_inventory.tsv` | Every HCD MS2 spectrum with all extracted metadata, its `precursor_status` and its `synthetic_record_id` |
| `flashlfq_input.tsv` | All synthetic FlashLFQ records (FlashLFQ "Generic" format) |
| `feature_mapping.tsv` | Spectrum → synthetic record → FlashLFQ feature, for every record FlashLFQ returned |
| `quantification_summary.json` | Parameters, image, per-file counts, effective FlashLFQ settings, status counts, QC statistics, and the pass/fail result of every validation check |
| `flashlfq/<run>/` | The per-run FlashLFQ input and the untouched FlashLFQ output (`QuantifiedPeaks.tsv`, `FlashLfqSettings.toml`, …) |
| `logs/pipeline.log`, `logs/extraction.log`, `logs/flashlfq_<run>.log` | Structured `key=value` logs. The FlashLFQ logs hold the full docker command and the container's stdout and stderr. |

### Main table columns

Required columns:

- `source_file`, `spectrum_title`, `native_scan_id`, `scan_number`
- `precursor_mz` (selected ion m/z), `precursor_charge`, `retention_time` (MS2 scan start time, minutes)
- `feature_id`, `feature_mz`, `feature_rt_start`, `feature_rt_apex`, `feature_rt_end`
- `ms1_intensity`, `quantification_status`

Additional columns:

- `internal_id` (`<file>::<nativeID>`), `spectrum_title_source`, `precursor_status`,
  `synthetic_record_id`, `neutral_mass`, `activation_method`, `precursor_intensity`,
  `thermo_monoisotopic_mz`
- `feature_charge` and `feature_charge_matches_precursor`
- `feature_n_charge_states`
- `feature_n_spectra`: the number of spectra that share this feature
- `feature_fwhm`, `feature_fwhm_status`, `feature_mass_error_ppm`
- `ms2_rt_minus_feature_apex` and `ms2_rt_within_feature_bounds`
- `flashlfq_detection_type`

`quantification_status` is one of:

| Status | Meaning | `ms1_intensity` |
|---|---|---|
| `quantified` | FlashLFQ found a chromatographic peak | FlashLFQ intensity (> 0) |
| `no_ms1_feature_found` | FlashLFQ returned the record without an apex: no envelope passed its mass, isotope and correlation filters. **This is not evidence that the precursor is absent.** | empty |
| `feature_with_zero_intensity` | The peak has an apex but zero intensity (not seen in practice) | 0 |
| `failed_flashlfq_run` | FlashLFQ failed or timed out for this run | empty |
| `failed_missing_from_flashlfq_output` | The record was sent but is missing from the FlashLFQ output (an integrity error) | empty |
| `not_quantified_<precursor_status>` | Not sent to FlashLFQ: `missing_precursor_charge`, `ambiguous_precursor_charge`, `missing_precursor_mz`, `ambiguous_multiple_precursors`, `ambiguous_multiple_selected_ions`, `missing_precursor`, `missing_selected_ion`, `invalid_precursor_charge`, `missing_or_invalid_retention_time` | empty |

An empty cell means missing. A measured zero is written as `0`.

### Joining with identifications

Join on **`source_file` + `spectrum_title`**. Titles are unique within a project (this is
checked). For search engines that report scan numbers instead of titles, join on the run name
(`source_file` without `.mzML[.gz]`) + `scan_number`. Several spectra can share one
`feature_id`. If you aggregate at the feature level later, deduplicate on `feature_id` first,
so that one intensity is not counted several times.

## How FlashLFQ is used (source analysis of FlashLFQ 2.1.5 / mzLib 1.0.592)

FlashLFQ has no identification-free mode. The CLI (`CMD/`, `Util/`) parses identification files.
The engine (`FlashLfqEngine`, in mzLib) quantifies `Identification` objects. The analysis below
shows that FlashLFQ can be used **unmodified**, by feeding it one synthetic Generic-format
record per MS2 event:

```
File Name  Scan Retention Time  Precursor Charge  Base Sequence  Full Sequence  Peptide Monoisotopic Mass  Protein Accession
<run>.mzML 19.468855            3                 NO_SEQUENCE    7a46706c58fa_SCAN1491  1802.760658         NO_PROTEIN_INFERENCE
```

FlashLFQ gives the two sequence columns different roles, and the pipeline fills them accordingly:

- **`Base Sequence` selects the isotope model, and nothing else that the pipeline uses.** It is the
  constant `NO_SEQUENCE` (point 2).
- **`Full Sequence` is the identity key.** FlashLFQ uses it for the isotope cache, grouping, peak
  comparison and the `QuantifiedPeaks.tsv` output. It holds the deterministic record id
  `<file hash>_SCAN<scan number>`. The file hash is the first 12 hex characters (lower case) of
  the SHA-1 of the run name, which is the file name without `.mzML` / `.mzML.gz`, UTF-8 encoded.
  Each file's hash is listed as `record_id_prefix` in `quantification_summary.json`. If a file
  has missing or repeated scan numbers, the mzML spectrum index is used instead
  (`<file hash>_INDEX5310`). The id is stored per spectrum as `synthetic_record_id` in all
  output tables. To regenerate it from a file name and a scan number:

  ```python
  hashlib.sha1("VEL180306-01_FWeiland180306H01_FW_180307072156".encode("utf-8")).hexdigest()[:12] + "_SCAN" + "1491"
  # '7a46706c58fa_SCAN1491'
  ```

1. **Mass.** FlashLFQ uses the supplied `Peptide Monoisotopic Mass` directly. The value is
   `z·(m/z − 1.007276466879)`, the same proton mass that mzLib uses.
2. **Isotope model.** `CalculateTheoreticalIsotopeDistributions` builds the formula from the
   `Base Sequence` only if every character is an amino-acid residue
   (`FlashLfqEngine.cs`, `SequenceResiduesAreValid(id.BaseSequence)`). `NO_SEQUENCE` contains `_`,
   so it is never valid. FlashLFQ therefore models each envelope with **averagine scaled to the
   supplied mass**, which is the neutral model for a precursor of unknown composition. Record
   ids made of amino-acid letters would change the quantification: on PXD002211, scan digits
   written as residues shifted the apex for 9% of spectra and biased intensities by up to ~7%
   (`validation/PXD002211_comparison/`).
3. **No unintended filtering.** The Generic reader filters on q-value and decoy status only
   when `Q-Value` or `Target/Decoy` columns exist, and the pipeline writes neither. The reader
   rejects a record if its sequence was already seen with a different mass, but every sequence
   is unique, so this never happens. No FASTA file, protein inference, normalisation or Bayesian
   statistics are used. `QuantifiedPeptides.tsv` and `QuantifiedProteins.tsv` are by-products
   and are ignored. The peptide table would also set ambiguous (merged) peaks to zero; the
   pipeline only reads the peak-level table.
4. **Uniqueness prevents unrelated merging.** Isotope distributions, peptide grouping and the
   mass-consistency check are all keyed by `Full Sequence`. With one unique ID per MS2 event,
   records are never grouped by sequence. They are merged only by FlashLFQ's own rule (below).
5. **One feature can serve several MS2 scans.** `RunErrorChecking` merges peaks that share the
   same apex MS1 peak (`MergeFeatureWith`). The merged peak keeps all its identifications.
6. **IDs survive into the output.** In `QuantifiedPeaks.tsv`, the `Full Sequence` column lists
   the distinct sequences of all merged identifications, joined with `|`, and `PSMs Mapped`
   gives their number (`458053a626eb_SCAN6332|458053a626eb_SCAN6445|458053a626eb_SCAN6994`).
   Record ids contain no `|`.
   The mapping is therefore by identifier,
   not by approximate m/z or RT. The pipeline checks the following:
   - every record appears in exactly one peak row
   - no ID is unknown or assigned to another run
   - `PSMs Mapped` equals the number of listed IDs
   - single-record rows echo the record's mass, charge and RT
   - the per-file row counts match the inventory

   Peaks without an apex are kept in the output as `no_ms1_feature_found`.
7. **Intensity definition.** An isotopic envelope's intensity is the sum of the observed
   isotope peaks, plus imputed values for expected isotopes that were not observed, **divided
   by the charge**. With `--integrate` (the pipeline default, FlashLFQ `--int`), the envelopes
   are summed over all MS1 scans of the peak, which approximates a peak area. With
   `--no-integrate` (FlashLFQ's own default), the intensity is the **apex envelope** only: the
   single MS1 scan with the highest envelope intensity. The validation confirmed the apex
   definition against the raw MS1 data.
8. **Charge states.** With `--id-specific-charge` (the pipeline default, FlashLFQ `--chg`),
   only the precursor's own charge state is searched. With `--no-id-specific-charge`
   (FlashLFQ's own default), every charge state between the minimum and maximum charge of the
   input is searched, and the reported apex can come from another charge state of the same
   neutral mass (12.8% of spectra in PXD009327; see `feature_charge_matches_precursor`). The
   same peptide fragmented at 2+, 3+ and 4+ then usually maps to one feature.

   Both defaults follow the comparison on PXD002211
   (`validation/PXD002211_comparison/PXD002211_v3_comparison_report.md`): with `--chg --int`,
   the results match the earlier ID-free quantification and the identification-based ionbot
   run, which use the same switches, most closely.
9. **Monoisotopic mass.** FlashLFQ compares each envelope with the ±1 Da shifted alternatives
   (`CheckIsotopicEnvelopeCorrelation`). If the selected m/z was not the monoisotopic peak, the
   envelope is usually rejected (`no_ms1_feature_found`) rather than quantified incorrectly. In
   PXD009327 the selected m/z always equals Thermo's monoisotopic m/z (ThermoRawFileParser).
   The trailer value is kept in `thermo_monoisotopic_mz` for other datasets.

### Two pitfalls in the FlashLFQ 2.1.5 CLI, and how the pipeline avoids them

- **Bool options are presence switches.** CommandLineParser treats `--nor false` as `--nor`
  plus an ignored token, so it **enables** normalisation. The same applies to `--int`, `--chg`,
  `--bay`, `--sil` and the others. The pipeline therefore never passes `true` or `false`; it
  passes a switch only when that switch must be on. It also checks the effective values in each
  run's `FlashLfqSettings.toml`, and a mismatch fails the run.
- **Match-between-runs cannot be turned off** (`--mbr` defaults to true). With stock defaults,
  FlashLFQ transferred synthetic records between runs: 103 MBR peaks in a 602-record test. The
  pipeline neutralises MBR with existing options, in two ways:
  - MBR donors need `QValue < pipfdr/5`. Synthetic records have QValue 0, so `--pipfdr 0`
    leaves **no donors**.
  - **Each run is quantified in its own FlashLFQ container**, so no other run exists to
    transfer from.

  The pipeline also checks that no `MBR` peak appears. Per-file runs and the joint run gave
  identical peak tables.

  One side effect remains: because MBR is enabled, FlashLFQ still writes a temporary peak index
  (`<run>.ind`) next to the spectra file. For this reason the spectra directory inside the
  container is a private, writable staging directory.

A one-line upstream fix would make `--mbr` usable (a nullable bool or `Default = false`). That
would require a custom image, and stock FlashLFQ already gives an equivalent result, so no
source change was made.

### Spectrum titles

The PXD009327 mzML files (ThermoRawFileParser 1.4.5) have **no** `MS:1000796` spectrum title.
When a title is present it is kept verbatim (`spectrum_title_source = mzml_spectrum_title`).
When it is absent, the pipeline builds one from `--title-template`, by default
`"<mzML run id>: <nativeID>"` (`spectrum_title_source = reconstructed`). For PXD009327 this
matches all 29408 `TITLE=` lines of the PRIDE MGF files exactly. If your PSMs come from MGF
files with a different title convention, set `--title-template`. The available placeholders are
`{run_id}`, `{native_id}`, `{file_stem}`, `{scan_number}` and `{index}`. Use
`--reference-mgf-dir` to check the result.

### HCD selection and precursor handling

A spectrum is eligible when all of the following are true:

- it has `ms level = 2`
- its activation has `MS:1000422` or `MS:1002481`
- it has no other dissociation term (CID, ETD, ECD, EThcD, …)

The `@hcd` filter string is used only when the mzML has no activation `cvParam` at all.
Non-HCD and MS3 spectra are counted in the summary but are not output rows. Spectra with
missing or ambiguous precursor information (several precursors or selected ions, no charge,
several possible charges, …) keep their row with an explicit status. The pipeline never picks
one candidate silently.

Retention times are converted to minutes from the `scan start time` unit (minute or second).
The original value and unit are kept in the inventory.

## Docker execution details

There is one `docker run` per mzML file:

- `--rm --network none --security-opt no-new-privileges`
- `--user $(id -u):$(id -g)`, so output files belong to you
- `--cpus <threads>` and an optional `--memory`
- `--name ms1quant-<pid>-<n>`, used by `docker kill` on timeout
- The spectra directory is mounted at `/data/mzml`. It is a private staging directory:
  - `.mzML.gz` files are decompressed there, because FlashLFQ reads only uncompressed mzML.
  - A plain `.mzML` input is bind-mounted **read-only** on top of it, so the original file
    cannot be modified.
- `output/flashlfq/<run>` is mounted at `/data/out`.
- `--mount` values are CSV-quoted, so paths containing spaces or commas work.

The FlashLFQ CLI exits with code 0 even when it aborts. The pipeline therefore also checks the
console for error messages, checks that `QuantifiedPeaks.tsv` exists, compares the number of
records read with the number sent, and checks the effective settings.

The image is pinned by digest. `smithchemwisc/flashlfq:2.1.5` is built from git commit
`754dbd98` (FlashLFQ tag 2.1.5). On an air-gapped server, use
`docker save` / `docker load` and keep the digest, or pass `--image` with a local tag.

## Validation

See `validation/PXD009327_validation_report.md` for the full report, including the server-side
checklist. In short, on the three PXD009327 runs:

- **29408 HCD spectra** (10238 / 9554 / 9616)
- 25049 quantified, 4359 `no_ms1_feature_found`, and 0 failures (current defaults, `--chg --int`)
- 2395 features shared by 5987 spectra
- record IDs regenerate from file name + scan number for all spectra
- all 18 pipeline checks passed
- apex peak found in the raw MS1 data for 300 of 300 sampled features
- titles match the PRIDE MGF files 1:1
- byte-identical results across repeated runs

FlashLFQ ran as the official 2.1.5 release build on the development machine, because no Docker
daemon was running there. Docker execution was then confirmed on a Linux server with
PXD002211 (6 runs, 163,789 spectra, all checks passed). See
`validation/PXD002211_comparison/PXD002211_v3_comparison_report.md` for the comparison with
the earlier ID-free quantification, and
`validation/PXD002211_comparison/PXD002211_record_id_comparison_report.md` for the comparison
with ionbot ID-based quantification, linked per scan by record id.

## Limitations

- Features are only as good as FlashLFQ's peak finding. For each feature, FlashLFQ requires
  ≥ 2 isotope peaks, a Pearson correlation > 0.7 with averagine, and the precursor's charge to be
  present in the peak. Low-abundance, high-charge or chimeric precursors are more often
  `no_ms1_feature_found`.
- Averagine is an approximation for non-peptide or heavily modified analytes, because the true
  composition is unknown when no identification is used.
- The intensity is the envelope intensity divided by z, summed over the peak's scans (`--int`, on
  by default), or the apex envelope only with `--no-integrate`. Intensities are not normalised
  between runs.
- The precursor mass comes from the mzML precursor m/z. When the Thermo "Monoisotopic M/Z"
  value is absent, this m/z is often not the monoisotopic peak. On PXD002211 this applies to
  49% of the scans, and 32% of the ionbot-identified PSMs have a precursor mass one or more
  isotopes off. These scans are quantified on a different envelope or not at all. The pipeline
  does not correct the monoisotopic precursor
  (`validation/PXD002211_comparison/PXD002211_record_id_comparison_report.md`).
- About 25% of quantified spectra have their MS2 RT outside FlashLFQ's peak bounds. This is
  mostly within one or two MS1 cycles, plus long-eluting background ions. The pipeline reports
  this and does not filter it.
- MBR is neutralised, not disabled (see above). FlashLFQ's MBR code path still runs, but it has
  no donors.
- Only uncompressed or gzip mzML is supported. Thermo RAW files would need FlashLFQ's licence
  handling and are not supported.

## Repository layout

```
quantify_precursors.py      pipeline (single file, standard library only)
requirements.txt            no runtime dependencies; pytest for tests
tests/                      pytest suite + fake_flashlfq.py (FlashLFQ stand-in)
validation/                 PXD009327 validation report, summary JSON, example output rows;
                            PXD002211 comparisons (PXD002211_comparison/) and server results
PXD009327-example/          one PXD009327 run (MZML, MGF, metadata) for the real-data test
```

The FlashLFQ source files cited in this README (`FlashLfqEngine.cs` and others) are in the
FlashLFQ 2.1.5 / mzLib source on GitHub. They are not part of this repository.

The hidden `--flashlfq-local-cmd "<dotnet> <CMD.dll>"` option runs a local FlashLFQ build
instead of Docker. It exists only for development and validation, and it logs a warning when
used.
