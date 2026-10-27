#!/usr/bin/env python
# coding: utf-8
import json, os, sys, re

HCD = {
    'accession': 'MS:1000422',
    'cvLabel': 'MS',
    'name': 'beam-type collision-induced dissociation',
    'value': 'HCD'
}

def json_HCD(path):
    try:
        with open(path) as user_file:
            return HCD in json.load(user_file)['ScanSettings']
    except:
        return False

def strip_compression(name):
    # remove .gz (current) or .xz (older projects) suffix
    return re.sub(r'\.(gz|xz)$', '', name)


for PXD in sys.argv[1:]:
    try:
        json_files = [_ for _ in os.scandir(os.path.join(PXD,'RAW'))]
    except:
            continue
    hcd_files  = []
    no_hcd_files = []
    for i,f in enumerate(json_files):
        PATH = os.path.join(PXD, 'metadata',
                            re.sub(r'\.raw$','-metadata.json',strip_compression(f.name),flags=re.I))
        if json_HCD(PATH):
            hcd_files.append(f)
        else:
            no_hcd_files.append(f)

    os.makedirs(os.path.join(PXD, 'metadata'), exist_ok=True)
    with open(os.path.join(PXD, 'metadata','CID_files.txt'), 'w') as outfile:
        # outfile.write(f'0\tCHECK\n')
        for i,f in enumerate(no_hcd_files):
            filename = strip_compression(f.name)
            outfile.write(f'{i+1}\t{filename}\n') 
    with open(os.path.join(PXD, 'metadata','HCD_files.txt'), 'w') as outfile:
        # outfile.write(f'0\tCHECK\n')
        for i,f in enumerate(hcd_files):
            filename = strip_compression(f.name)
            outfile.write(f'{i+1}\t{filename}\n') 
    print(os.path.join(PXD, 'metadata','HCD_files.txt'), len(hcd_files))
    print(f'{len(hcd_files)} HCD files identified')
    print(f'{len(no_hcd_files)} other files identified')
