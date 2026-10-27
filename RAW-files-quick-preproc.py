#!/usr/bin/env python
# coding: utf-8
import os, sys, argparse, subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

parser = argparse.ArgumentParser()
parser.add_argument('PXDs', nargs='+')
parser.add_argument('-j', '--jobs', type=int, default=4,
                    help='number of RAW files processed in parallel (default: 4)')
args = parser.parse_args()

def do_command(CMD):
    # no pipe: output goes straight to the terminal (or log file), as it is written
    print(">",CMD, flush=True)
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

to_process = args.PXDs
j = len(to_process)
for i,PXD in enumerate(to_process):
    print('\n',PXD,'\t',f'{i+1}/{j}')

    raw_dir  = os.path.join(PXD,'RAW')
    meta_dir = os.path.join(PXD,'metadata')
    mgf_dir  = os.path.join(PXD,'MGF')
    mzml_dir = os.path.join(PXD,'mzml')
    if not os.path.isdir(raw_dir):
        print("WARNING! -- no RAW files here!")
        do_command(f"touch {os.path.join(PXD,'no-hcd-files')}")
        continue

    if not os.path.isdir(meta_dir):
        decompress_raw(raw_dir)
        do_command(CMD = f'bash get-rawfile-metadata.sh {PXD} {args.jobs}')
    if os.path.isdir(meta_dir):
        do_command(CMD = f'python parse-json-metadata.py {PXD}')

    PATH = os.path.join(PXD,'metadata','HCD_files.txt')
    try:
        hcd_files = list(pd.read_csv(PATH, sep='\t', names=['Files'], index_col=0).Files)
    except:
        hcd_files = []
    if len(hcd_files) < 1:
        print(f'Warning! -- No HCD files in {PXD}!')
        do_command(f"touch {os.path.join(PXD,'no-hcd-files')}")


    # do_command(CMD = f'gzip {raw_dir}/*.raw {raw_dir}/*.RAW')