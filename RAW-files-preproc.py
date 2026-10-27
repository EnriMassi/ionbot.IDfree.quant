#!/usr/bin/env python
# coding: utf-8
import os, sys, argparse, subprocess
from concurrent.futures import ThreadPoolExecutor

parser = argparse.ArgumentParser()
parser.add_argument('PXDs', nargs='+')
parser.add_argument('-j', '--jobs', type=int, default=4,
                    help='number of RAW files processed in parallel (default: 4)')
args = parser.parse_args()

def do_command(CMD):
    # no pipe: output goes straight to the terminal (or log file), as it is written
    print(">", CMD, flush=True)
    subprocess.run(CMD, shell=True)

def run_parallel(CMDs):
    with ThreadPoolExecutor(args.jobs) as executor:
        list(executor.map(do_command, CMDs))

def decompress_raw(raw_dir):
    # RAW files are now stored as .gz; .xz is kept for older projects
    CMDs = []
    for f in os.scandir(raw_dir):
        if os.path.exists(f.path[:-3]): continue
        if f.name.endswith('.gz'):
            CMDs.append(f'gzip -d {f.path}')
        elif f.name.endswith('.xz'):
            CMDs.append(f'unxz -T16 {f.path}')
    run_parallel(CMDs)

def compress_raw(raw_dir):
    # only compress files that are not compressed yet (.gz, or legacy .xz)
    # -f overwrites a partial .gz left behind by an interrupted run
    run_parallel([f'gzip -f {f.path}' for f in os.scandir(raw_dir)
                  if f.is_file() and not f.name.endswith(('.gz','.xz'))])

for PXD in args.PXDs:
    raw_dir  = os.path.join(PXD,'RAW')
    meta_dir = os.path.join(PXD,'metadata')
    mgf_dir  = os.path.join(PXD,'MGF')
    mzml_dir = os.path.join(PXD,'mzml')

    if os.path.isfile(os.path.join(PXD,'flashlfq.tsv')):
        if os.path.isdir(raw_dir):
            compress_raw(raw_dir)
        continue

    if not os.path.isdir(raw_dir):
        print("WARNING! -- no RAW files here!")
        do_command(f"touch {os.path.join(PXD,'no-hcd-files')}")
        continue

    if not os.path.isdir(meta_dir):
        decompress_raw(raw_dir)
        do_command(CMD = f'bash get-rawfile-metadata.sh {PXD} {args.jobs}')
    if os.path.isdir(meta_dir):
        do_command(CMD = f'python parse-json-metadata.py {PXD}')
    # get-mgf.sh and get-mzml.sh both decompress RAW files: do not run them at the same time
    if not os.path.isdir(mgf_dir):
        do_command(CMD = f'bash get-mgf.sh {PXD} {args.jobs}')
    if not os.path.isdir(mzml_dir):
        do_command(CMD = f'bash get-mzml.sh {PXD} {args.jobs}')

    compress_raw(raw_dir)
