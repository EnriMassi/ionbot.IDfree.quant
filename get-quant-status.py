#!/usr/bin/env python
# coding: utf-8
import json, os, sys, re
import pandas as pd
import numpy as np

projects = pd.read_csv('pxds.to.quantify.txt', sep='\t', names=['path'])
print(projects.head())
print(projects.tail())

def is_quant_successful(project_path):
    return int(os.path.exists(os.path.join(project_path,'QuantifiedPeaks.tsv')) or os.path.exists(os.path.join(project_path,'QuantifiedPeaks_v2.tsv')))

def is_quant_done(project_path):
    return int(os.path.exists(os.path.join(project_path,'QuantifiedPeaks.tsv')) or 
               os.path.exists(os.path.join(project_path,'QuantifiedPeaks_v2.tsv')) or 
               os.path.exists(os.path.join(project_path,'ErrorReport.txt')) or 
               os.path.exists(os.path.join(project_path,'ErrorReport_v2.txt')))

def is_quant_done_v2(project_path):
    return int(os.path.exists(os.path.join(project_path,'QuantifiedPeaks_v2.tsv')) or os.path.exists(os.path.join(project_path,'ErrorReport_v2.txt')))

def has_NO_hcd_files(project_path):
    return int(os.path.exists(os.path.join(project_path,'no-hcd-files')))

def count_hcd_files(project_path):
    if has_NO_hcd_files(project_path):
        return 0 
    try:
        hcd_files = pd.read_csv(os.path.join(project_path,'metadata','HCD_files.txt'), sep='\t', names=['Files'], index_col=0)
        return len(hcd_files)
    except:
        return np.nan

def count_meta_files(project_path):
    try:
        meta_files = [_ for _ in os.scandir(os.path.join(project_path,'metadata')) if _.name.endswith('.json')]
        return len(meta_files)
    except:
        return 0
    
def count_files(project_path):
    try:
        files = [_ for _ in os.scandir(os.path.join(project_path,'RAW'))]
        return len(files)
    except:
        return 0


projects['N_raw_files']   = projects.path.apply(count_files)
projects['N_meta_files']  = projects.path.apply(count_meta_files)
projects['N_hcd_files']   = projects.path.apply(count_hcd_files)
projects['is_quant_done'] = projects.path.apply(is_quant_done)
projects['successful']    = projects.path.apply(is_quant_successful)  
projects['is_quant_done_v2']    = projects.path.apply(is_quant_done_v2) 

projects.sort_values('N_hcd_files', inplace=True)

projects.to_csv('Quantification-status.csv', index=False)
