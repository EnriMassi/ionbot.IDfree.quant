# ionbot.IDfree.quant

Identification-free MS1 quantification of Thermo RAW files from PRIDE projects.

For every HCD spectrum in a project, the pipeline measures the MS1 intensity of its precursor with [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ). No database search is needed: each MS2 spectrum from the MGF file is given to FlashLFQ as a placeholder "peptide" with its own retention time, charge and precursor mass. The quantified peaks can later be linked to identifications through the spectrum title.

## Requirements

- Linux with bash, `gzip` and `xz` (`unxz` is only needed for projects compressed with xz in the past).
- Docker, with these images:
  - `mono`, to run ThermoRawFileParser: `docker pull mono`
  - `smithchemwisc/flashlfq:1.2.4`: `docker pull smithchemwisc/flashlfq:1.2.4`
- [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) 1.3.4, unpacked into a folder called `ThermoRawFileParser/` inside this repository.
- Python 3 with `pandas`, `numpy` and `pyteomics`. `combine-quant-results.py` also needs `polars`; `walk-through-pride-archive.py` needs `requests` and `beautifulsoup4`.

The docker commands mount `/public/` into the containers and use the current directory as working directory. Therefore this repository and all project folders must be located under `/public/`. Change the `-v /public/:/public/` option in the scripts if your data are elsewhere.

All scripts must be started from the root of this repository, because they call `ThermoRawFileParser/` and `ionbot.quant/` with relative paths.

## Input

One folder per PRIDE project, with the RAW files in a `RAW` subfolder. The RAW files can be compressed with gzip (`.raw.gz`) or xz (`.raw.xz`), or uncompressed:

```
PXD000001/
└── RAW/
    ├── sample_01.raw.gz
    └── sample_02.raw.gz
```

## Usage

List the project folders to process in `pxds.to.quantify.txt`, one path per line, then run:

```bash
bash quant-from-mgfs.sh [n_pxds] [n_jobs]
```

- `n_pxds` (default 2): number of projects processed at the same time.
- `n_jobs` (default 4): number of RAW files converted in parallel within each project.

Up to `n_pxds * n_jobs` ThermoRawFileParser containers can run at once. FlashLFQ always runs on one file at a time within each project. The output of each project goes to `<PXD>/quant.log`.

Example, with 3 projects at a time and 4 RAW files per project:

```bash
bash quant-from-mgfs.sh 3 4
```

## Pipeline steps

`quant-from-mgfs.sh` runs these two steps for every project. Both steps can also be run on their own.

### 1. Preprocessing: `RAW-files-preproc.py`

```bash
python RAW-files-preproc.py [-j N_JOBS] PXD [PXD ...]
```

1. Decompresses the RAW files.
2. Extracts the metadata of each RAW file as JSON (`get-rawfile-metadata.sh`).
3. Sorts the RAW files into HCD and non-HCD files from the metadata (`parse-json-metadata.py`). The lists are written to `metadata/HCD_files.txt` and `metadata/CID_files.txt`.
4. Converts the HCD files to MGF (`get-mgf.sh`) and to mzML (`get-mzml.sh`).
5. Compresses the RAW files again with gzip.

Each step is skipped if its output folder already exists, so an interrupted project can be restarted with the same command. The conversion scripts write into a temporary folder (`metadata.tmp`, `MGF.tmp`, `mzml.tmp`) and rename it only when all files are done. After an interruption, only the temporary folder is left behind, and the next run rebuilds it.

`RAW-files-quick-preproc.py` runs only steps 1 to 3. Use it to find out which projects contain HCD files before converting them.

### 2. Quantification: `ionbot.quant/mgf2FlashLFQ.py`

```bash
python ionbot.quant/mgf2FlashLFQ.py PXD [PXD ...]
```

This script is adapted from `ionbot2FlashLFQ.py` in [sdgroeve/ionbot.quant](https://github.com/sdgroeve/ionbot.quant) (https://github.com/sdgroeve/ionbot.quant). The original script builds the FlashLFQ input file from ionbot identifications; this version builds it from the MGF spectra instead. Refer to the original repository for the FlashLFQ input format and docker usage when editing this step.

For every HCD file of a project:

1. Reads the MGF file and writes a FlashLFQ input file with one entry per spectrum. The precursor mass is computed as `(m/z - 1.007276) * charge`. The spectrum title is stored in the `Protein Accession` column.
2. Runs FlashLFQ on the corresponding mzML file (`--chg --int --ppm=20`).

The results of all files are combined into `<PXD>/QuantifiedPeaks.tsv`. The `Protein Group` column of this file contains the spectrum titles. Projects that already have `QuantifiedPeaks.tsv` or `ErrorReport.txt` are skipped.

## Output

```
PXD000001/
├── RAW/                      RAW files, compressed again with gzip
├── metadata/                 JSON metadata, HCD_files.txt, CID_files.txt
├── MGF/                      MGF files (.mgf.gz)
├── mzml/                     mzML files (.mzML.gz)
├── QuantifiedPeaks.tsv       FlashLFQ peaks of all HCD files
└── quant.log                 log of quant-from-mgfs.sh
```

Other files that can appear in a project folder:

- `no-hcd-files`: the project has no RAW files, or none of them is HCD.
- `ErrorReport.txt`: FlashLFQ produced no results for the project.

## Utility scripts

| Script | Purpose |
|---|---|
| `walk-through-pride-archive.py` | Crawls the PRIDE FTP archive and writes the list of file URLs to `PRIDE-URLs-files.csv`. |
| `Download-from-PRIDE.sh` | Downloads the `.raw` files of one PRIDE project: `bash Download-from-PRIDE.sh <project FTP URL> <target folder>`. |
| `decompress-raw.sh` | Decompresses the HCD RAW files of one project: `bash decompress-raw.sh <PXD>`. |
| `double-check-metadata.sh` | Extracts the metadata of the files listed in `metadata/CID_files.txt` again, to check files that were not classified as HCD. |
| `ionbot.quant/ionbot2FlashLFQ-project.py` | Quantifies the ionbot identifications of a project the same way as `mgf2FlashLFQ.py`, for comparison: `python ionbot.quant/ionbot2FlashLFQ-project.py [-x] [-i IONBOT_DIR] PXD [PXD ...]`. Reads `<PXD>/IONBOT_v*/<file>.mgf.gzip/ionbot.first.csv` (highest version, or the folder given with `-i`) and the mzML files from `mzml/`, and writes `<PXD>/QuantifiedPeaks.ionbot.tsv` (or `ErrorReport.ionbot.txt`). |
| `get-quant-status.py` | Writes `Quantification-status.csv` with the number of RAW, metadata and HCD files and the quantification status of each project in `pxds.to.quantify.txt`. |
| `ionbot.quant/combine-quant-results.py` | Concatenates the TSV files listed in a text file into one CSV: `python ionbot.quant/combine-quant-results.py paths.txt -o combined.csv`. |

## Stopping a run

Press Ctrl+C or send `kill` to the process. The conversion scripts stop their running docker containers before they exit. Check with `docker ps` that no ThermoRawFileParser or FlashLFQ container is left, and stop any leftover container with `docker stop <container id>`.

## Acknowledgements

The FlashLFQ step in `ionbot.quant/` is based on [sdgroeve/ionbot.quant](https://github.com/sdgroeve/ionbot.quant), which converts ionbot search results into FlashLFQ input files.

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE).
