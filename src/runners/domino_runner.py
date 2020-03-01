#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys
sys.path.insert(0, '../')

import os
import numpy as np

from src.utils import infra
from src import constants
from src.implementations.domino import main as domino_main
from src.utils.network import build_all_reports
import src.utils.server as server
from src.utils.ensembl2entrez import ensembl2entrez_convertor
ALGO_NAME = "domino_original"

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

def init_specific_params(score_file_name, dest_algo_dir):

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)

    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    last_q_index = np.where(deg_data[:,np.where(h_cols=="qval")[0][0]]>0.05)[0][0]
    ge_list_file_name = os.path.join(constants.OUTPUT_DIR, "ge_list.txt")
    open(ge_list_file_name, "w+").write("\n".join([x for x in h_rows[:last_q_index] if len(ensembl2entrez_convertor([x]))>0 ]))
    return ge_list_file_name


def main(active_genes_file, network_file_name, output_folder, score_method=constants.PREDEFINED_SCORE):
    print("run domino....")

    network_file_name, score_file_name, bg_genes = server.init_common_params(network_file_name, active_genes_file=active_genes_file)

    modules=domino_main(active_genes_file=active_genes_file, network_file=network_file_name, slices_file=os.path.join(constants.NETWORKS_DIR, "dip_ng_modularity_components.txt"), slice_threshold=module_sig_th)

    modules = filter(lambda x: len(x) > 3, modules)
    all_bg_genes = [bg_genes for x in modules]

    if constants.REPORTS:
        build_all_reports(ALGO_NAME, output_folder, modules, all_bg_genes, score_file_name, network_file_name)


if __name__ == "__main__":
    ds=["PASCAL_SUM_Crohns_Disease.G50"] # ["PASCAL_SUM_Breast_Cancer.G50", "PASCAL_SUM_Crohns_Disease.G50", "PASCAL_SUM_Schizophrenia.G50", "PASCAL_SUM_Triglycerides.G50", "PASCAL_SUM_Type_2_Diabetes.G50"]
    for cur in ds:
        constants.update_dirs(DATASET_NAME_u=cur) # Type_2_Diabetes Crohns_Disease
        main(dataset_name=constants.DATASET_NAME, score_method=constants.PREDEFINED_SCORE, module_sig_th=0.3)



