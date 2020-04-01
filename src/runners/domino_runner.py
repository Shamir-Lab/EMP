#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys
sys.path.insert(0, '../')

import os
import pandas as pd

from src import constants
from src.implementations.domino import main as domino_main
from src.utils.network import build_all_reports
from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.network import get_network_genes

ALGO_NAME = "DOMINO"

def extract_modules_and_bg(bg_genes, dest_algo_dir):
    results = open(os.path.join(dest_algo_dir, "modules.txt")).readlines()
    modules = [[] for x in range(max([int(x.strip().split(" =")[1]) for x in results[1:]]) + 1)]
    for x in results[1:]:
        if int(x.strip().split(" =")[1]) != -1:
            modules[int(x.strip().split(" =")[1])].append(x.strip().split(" =")[0])
        else:
            modules.append([x.strip().split(" =")[0]])
    modules = filter(lambda x: len(x) > 3, modules)
    all_bg_genes = [bg_genes for x in modules]
    print("extracted {} modules".format(len(modules)))
    return modules, all_bg_genes

def init_params(dataset_file_name, network_file_name, output_folder):


    df_scores=pd.read_csv(dataset_file_name, sep='\t', index_col=0)
    sig_genes=df_scores['qval'][df_scores['qval']<0.05].index
    active_genes_file=os.path.join(output_folder, "active_genes_file.txt")
    open(active_genes_file, "w+").write("\n".join([x for x in sig_genes if len(ensembl2entrez_convertor([x]))>0 ]))
    bg_genes=get_network_genes(network_file_name)
    return active_genes_file, bg_genes


def main(dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):
    print("run domino....")
    slices_file=kwargs['slices_file']

    slice_threshold=0.3
    if 'slice_threshold' in kwargs:
        slice_threshold=kwargs['slice_threshold']

    module_threshold=0.05
    if 'module_threshold' in kwargs:
        module_threshold=kwargs['module_threshold']

    active_genes_file, bg_genes = init_params(dataset_file_name, network_file_name, output_folder)

    modules=domino_main(active_genes_file=active_genes_file, network_file=network_file_name, slices_file=slices_file, slice_threshold=slice_threshold, module_threshold=module_threshold)

    modules = list(filter(lambda x: len(x) > 3, modules))
    all_bg_genes = [bg_genes for x in modules]

    build_all_reports(ALGO_NAME, modules, all_bg_genes, go_folder, os.path.join(output_folder, "report"))


if __name__ == "__main__":
    ds=["PASCAL_SUM_Crohns_Disease.G50"] # ["PASCAL_SUM_Breast_Cancer.G50", "PASCAL_SUM_Crohns_Disease.G50", "PASCAL_SUM_Schizophrenia.G50", "PASCAL_SUM_Triglycerides.G50", "PASCAL_SUM_Type_2_Diabetes.G50"]
    for cur in ds:
        constants.update_dirs(DATASET_NAME_u=cur) # Type_2_Diabetes Crohns_Disease
        main(dataset_name=constants.DATASET_NAME, score_method=constants.PREDEFINED_SCORE, module_sig_th=0.3)



