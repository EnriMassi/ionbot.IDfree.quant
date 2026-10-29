#!/usr/bin/env python
# coding: utf-8
import sys, os, re, shutil, subprocess, hashlib
import pandas as pd
import numpy as np
from pyteomics import mgf

PROTON_MASS = 1.007276

def remove_extension(x):
    return ".".join(x.split(".")[:-1])

def do_command(CMD):
    # no pipe: output goes straight to the terminal (or log file), as it is written
    print(">",CMD, flush=True)
    subprocess.run(CMD, shell=True)

def decompress(path):
    # path is the uncompressed file; .gz is current,
    # .gzip (MGF) and .xz (mzML) come from older runs
    if os.path.exists(path): return
    if os.path.exists(path+'.gz'):
        do_command(f'gzip -d {path}.gz')
    elif os.path.exists(path+'.gzip'):
        do_command(f'gzip -S .gzip -d {path}.gzip')
    elif os.path.exists(path+'.xz'):
        do_command(f'unxz -T16 {path}.xz')

def compress(path):
    # -f overwrites a partial .gz left behind by an interrupted run
    if os.path.exists(path):
        do_command(f'gzip -f {path}')

def digit_to_aminoacid(x):
    aa_symbols = 'ACDEFGHILM'
    x = str(x)
    aa = [aa_symbols[int(_)] for _ in x]
    return ''.join(aa)

for PXD in sys.argv[1:]:
    # skip if we already have quant
    if os.path.exists(os.path.join(PXD,'QuantifiedPeaks_v4.tsv')) or os.path.exists(os.path.join(PXD,'ErrorReport_v4.txt')):
        continue

    # unzip stuff
    PATH = os.path.join(PXD,'metadata','HCD_files.txt')
    try:
        hcd_files = list(pd.read_csv(PATH, sep='\t', names=['Files'], index_col=0).Files)
    except:
        hcd_files = []
    if len(hcd_files) < 1:
        print(f'Warning! -- No HCD files in {PXD}!')
        do_command(f"touch {os.path.join(PXD,'no-hcd-files')}")
        print(PXD,'end')
        continue
    else:
        print(f'# HCD files = {len(hcd_files)}')

    # FlashLFQ is given the whole mzml folder: compress the mzML files left
    # by RAW-files-preproc.py, so that only the current one is uncompressed
    mzml_dir = os.path.join(PXD,'mzml')
    if os.path.isdir(mzml_dir):
        for _ in os.scandir(mzml_dir):
            if _.name.endswith('.mzML'):
                compress(_.path)

    # read mgfs
    for i,f in enumerate(hcd_files):
        print(i+1,'/',len(hcd_files))
        raw_name  = re.sub(r'\.raw$', '', f, flags=re.I)
        mgf_path  = os.path.join(PXD, 'MGF',  raw_name+'.mgf')
        mzml_path = os.path.join(PXD, 'mzml', raw_name+'.mzML')

        fname  = remove_extension(os.path.split(mgf_path)[-1])
        outdir = os.path.join(PXD,'FlashLFQ',fname)

        if not os.path.exists(os.path.join(outdir,'QuantifiedPeaks.tsv')):
            os.makedirs(outdir, exist_ok=True)
            flashlfqpath = os.path.join(outdir,'flashlfq.tsv')
            decompress(mgf_path)
            decompress(mzml_path)
            tmp = []
            
            try:
                mgf.read(mgf_path)
            except:
                compress(mgf_path)
                compress(mzml_path)
                print("------------")
                continue

            for _ in mgf.read(mgf_path):
                try:
                    rt     = _['params']['rtinseconds'] / 60
                    charge = _['params']['charge'][0]
                    mz     = _['params']['pepmass'][0]
                    title  = _['params']['title']

                    scan = int(re.search(r'scan=(\d+)', title).group(1))

                    # tmp.append([rt, charge, 
                    #             (mz-PROTON_MASS)*charge, 
                    #             title, fname, False, scan_n])
                    tmp.append(
                        {
                            'Scan Retention Time': rt,
                            'Precursor Charge': charge,
                            'Peptide Monoisotopic Mass': (mz-PROTON_MASS)*charge,
                            'Protein Accession': title,
                            'File Name': fname,
                            'Shared_pep': False,
                            'Base Sequence': digit_to_aminoacid(scan),
                            'Full Sequence': hashlib.sha1(fname.encode("utf-8")).hexdigest()[:12] + "_SCAN" + str(scan)
                        }
                    )
                except:
                    continue
            tmp = pd.DataFrame(tmp)
            tmp.to_csv(flashlfqpath, sep="\t", index=False)

            try:
                do_command(f"docker run --rm -v /public/:/public/ -w {os.getcwd()} smithchemwisc/flashlfq:1.2.4 "+
                      f"--idt {flashlfqpath} --rep {PXD}/mzml --chg --int --ppm=20 --out={outdir}")
            except:
                print('problem with files(?)')
        compress(mgf_path)
        compress(mzml_path)
        print("------------")

    combo_quant_peaks = []
    for _ in os.scandir(os.path.join(PXD,'FlashLFQ')):
        try:
            combo_quant_peaks.append(pd.read_csv(os.path.join(_.path,'QuantifiedPeaks.tsv'), sep='\t', low_memory=False))
        except:
            pass
    print(f'Combining {len(combo_quant_peaks)} output files...')
    try:
        pd.concat(combo_quant_peaks, ignore_index=True).to_csv(os.path.join(PXD,'QuantifiedPeaks_v4.tsv'), sep='\t', index=False)
    except:
        do_command(f"touch {os.path.join(PXD,'ErrorReport_v4.txt')}")
    shutil.rmtree(os.path.join(PXD,'FlashLFQ'))
    print("------------")