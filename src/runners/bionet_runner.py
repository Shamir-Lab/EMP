#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../..')

import os
import time
import shutil
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()



import src.constants as constants

from src.utils.r_runner import run_rscript
from src.utils.network import remove_subgraph_by_nodes
from src.utils.network import build_all_reports
import numpy as np


ALGO_NAME = "bionet"
ALGO_DIR = "/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/bionet"


def run_bionet(deg_file_name, network_file_name, fdr=0.05):
    script = open(os.path.join(ALGO_DIR, "{}.r".format(ALGO_NAME))).read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name, fdr=fdr)


def init_specific_params(network_file_name, omitted_genes = [], ts=str(time.time())):
    return remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)


def get_module(network_file_name, score_file_name, omitted_genes, ts=str(time.time()),fdr=0.05):
    network_file_name = init_specific_params(network_file_name=network_file_name, omitted_genes=omitted_genes)
    results = run_bionet(score_file_name, network_file_name,fdr)
    module_genes = np.array(results["module_genes"])

    bg_genes = results["bg_genes"]
    # open(os.path.join(constants.OUTPUT_DIR, "{}_module_genes_{}.txt".format(ALGO_NAME, ts)), "w+").write(
    #     "\n".join(module_genes))
    # open(os.path.join(constants.OUTPUT_DIR, "{}_bg_genes_{}.txt".format(ALGO_NAME, ts)), "w+").write(
    #     "\n".join(bg_genes))

    sys.stdout.write("module gene size: {}. ratio: {}\n".format(module_genes.shape[0], module_genes.shape[0]/float(bg_genes.shape[0])))
    return list(module_genes), list(bg_genes)


def run_bionet_for_all_modules(fdr, network_file_name, score_file_name):
    omitted_genes = []
    modules = []
    all_bg_genes = []
    small_modules=0
    for x in range(50):
        modules_genes=None
        bg_genes=None
        try:
           module_genes, bg_genes = get_module(network_file_name, score_file_name, omitted_genes, str(x), fdr=fdr)
        except Exception as e:
           print("got an exception while trying to extract module #{}:\n{}".format(x, e))
           break
        if len(module_genes) == 0 or small_modules==5: break

        omitted_genes += list(module_genes)
        if len(module_genes) > 3:
            small_modules=0
            modules.append(module_genes)
            all_bg_genes.append(bg_genes)
        else:
            small_modules += 1
    return all_bg_genes, modules

def main(dataset_file_name, network_file_name, go_folder, output_folder, fdr=0.05, **kwargs):

    all_bg_genes, modules = run_bionet_for_all_modules(fdr, network_file_name, dataset_file_name)

    build_all_reports(ALGO_NAME, modules, all_bg_genes, go_folder, os.path.join(output_folder, "report"))




if __name__ == "__main__":

    main(dataset_file_name="/specific/netapp5/gaga/hagailevi/emp_test/original_datasets/tnfa.tsv", network_file_name="/specific/netapp5/gaga/hagailevi/emp_test/networks/dip.sif", go_folder="/specific/netapp5/gaga/hagailevi/emp_test/networks/dip.sif", output_folder="/specific/netapp5/gaga/hagailevi/emp_test/true_solutions/")







