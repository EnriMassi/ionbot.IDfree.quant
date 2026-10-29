#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
import sys
import polars as pl

"""
combine-quant-results.py

Read a text file listing paths to .tsv files (one per line), read them with polars,
concatenate into a single DataFrame, and write out as CSV.
"""


def read_paths(list_file: Path):
    with list_file.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            yield Path(line)

def parse_cli():
    p = argparse.ArgumentParser(description="Concatenate TSV files listed in a file into one CSV")
    p.add_argument("list_file", type=Path, help="Text file with paths to .tsv files (one per line)")
    p.add_argument("-o", "--out", type=Path, default=Path("combined.csv"), help="Output CSV path")
    p.add_argument("--sep", default="\t", help="Input TSV separator (default: tab)")
    p.add_argument("--skip-missing", action="store_true", help="Skip missing files instead of failing")
    return p.parse_args()

def main():
    args = parse_cli()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not args.list_file.exists():
        logging.error("List file not found: %s", args.list_file)
        sys.exit(2)

    frames = []
    my_files = list(read_paths(args.list_file))
    total_files = len(my_files)
    logging.info("Processing %d files from %s", total_files, args.list_file)
    for idx,fp in enumerate(my_files):
        if not fp.exists():
            msg = f"File not found: {fp}"
            if args.skip_missing:
                logging.warning(msg + " (skipped)")
                continue
            else:
                logging.error(msg)
                sys.exit(3)
        try:
            logging.info(f"Reading {fp} : file {idx + 1} of {total_files}")
            df = pl.read_csv(
                str(fp), 
                separator=args.sep, 
                null_values=["NA", "NaN", "", '-'],
                schema_overrides={"File Name": pl.String}
                )
            frames.append(df)
        except Exception as e:
            logging.error("Failed to read %s: %s", fp, e)
            sys.exit(4)

    if not frames:
        logging.error("No data frames were loaded. Exiting.")
        sys.exit(5)

    logging.info("Concatenating %d frames", len(frames))
    out_df = pl.concat(frames, how="vertical", rechunk=True)

    logging.info("Writing combined CSV to %s", args.out)
    out_df.write_csv(str(args.out))
    logging.info("Done")


if __name__ == "__main__":
    main()