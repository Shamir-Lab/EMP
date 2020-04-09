#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../..')

import os
import numpy as np
import subprocess

import shutil
import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import build_all_reports
from src.utils.network import remove_subgraph_by_nodes

from src.utils.network import get_network_genes

import src.utils.infra as infra

from src.utils.network import output_modules

ALGO_NAME = "keypathwayminer"
ALGO_DIR = "/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/keypathwayminer"

NETWORK_NAME = "dip"

import random

def init_params(score_file_name, omitted_genes, network_file_name, ts, dest_algo_dir):
    if os.path.exists(os.path.join(dest_algo_dir, "results")):
        shutil.rmtree(os.path.join(dest_algo_dir, "results"))

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    sig_binary_col = deg_data[:,np.where(h_cols=="qval")[0][0]]<0.05
    sig_binary_output = np.c_[h_rows,  np.array(sig_binary_col, dtype=np.int)]
    score_file_name = os.path.join(dest_algo_dir, "binary_score.txt")
    open(score_file_name, "w+").write("\n".join(["\t".join(["id", "pval", "qval"])] + ["\t".join(list(x) + list([x[-1]])) for x in sig_binary_output]))


    new_network_file_name = remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)

    bg_genes=get_network_genes(network_file_name)

    return score_file_name, new_network_file_name, bg_genes


def format_scripts(score_file_name, network_name="dip", STRATEGY="INES", algorithm="GREEDY", algo_dir=None, dataset_name=None):
    script_file_name=format_script(os.path.join(constants.dir_path,"src","sh", "run_{}.sh".format(ALGO_NAME)), STRATEGY=STRATEGY, ALGORITHM=algorithm)
    format_script(os.path.join(algo_dir, "kpm.properties"), uid=False, network_name=network_name, algorithm=algorithm)
    format_script(os.path.join(algo_dir, "datasets_file.txt"), uid=False, score_file_name=score_file_name)

    return script_file_name


def extract_module_genes(bg_genes, STRATEGY, algorithm, dest_algo_dir):
    i = 1
    modules = []
    while os.path.exists(os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))) and i<2:
        results = open(
            os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))).readlines()
        results = list(map(lambda x: x.strip(), results))
        modules.append(results)
        i += 1

    # module_genes = list(set([y for x in modules for y in x]))
    # module_genes = list(set(module_genes))
    # module_genes_file_name = os.path.join(constants.OUTPUT_DIR, "{}_{}_{}_module_genes.txt".format(ALGO_NAME, STRATEGY, algorithm))
    # open(module_genes_file_name, "w+").write("\n".join(module_genes))

    return modules, [bg_genes for x in modules]


def main(dataset_file_name, network_file_name, go_folder, output_folder, fdr=0.05, **kwargs):
    strategy = "INES"
    algorithm = "GREEDY"
    omitted_genes = []
    modules = []
    all_bg_genes = []
    dest_algo_dir = "{}_{}".format(ALGO_DIR, random.random())
    shutil.copytree(ALGO_DIR, dest_algo_dir)
    empty_counter = 0
    for cur_i_module in range(50):
        binary_score_file_name, cur_network_file_name, bg_genes= init_params(dataset_file_name, omitted_genes,
                                                                         network_file_name, str(random.random()), dest_algo_dir)

        script_file_name=format_scripts(score_file_name=binary_score_file_name, network_name=cur_network_file_name,
                       STRATEGY=strategy, algorithm=algorithm, algo_dir=dest_algo_dir)
        print(subprocess.Popen("bash {}".format(script_file_name), shell=True,
                               stdout=subprocess.PIPE, cwd=dest_algo_dir).stdout.read())
        module, all_bg_gene = extract_module_genes(bg_genes, strategy, algorithm, dest_algo_dir)
        print(module)
        if len(module[0]) > 3:
            empty_counter=0
            modules.append(module[0])
            all_bg_genes.append(all_bg_gene[0])
        else:
            empty_counter+=1
        omitted_genes += list(module[0])
        os.remove(script_file_name)

        if empty_counter>3:
            print("got more that 3 smalle modules in row. continue...")
            break

    shutil.rmtree(dest_algo_dir)

    build_all_reports(ALGO_NAME, modules, all_bg_genes, go_folder, os.path.join(output_folder, "report"))



if __name__ == "__main__":

    main(dataset_file_name="/specific/netapp5/gaga/hagailevi/emp_test/original_datasets/brca.tsv", network_file_name="/specific/netapp5/gaga/hagailevi/emp_test/networks/dip.sif", go_folder="/specific/netapp5/gaga/hagailevi/emp_test/networks/dip.sif", output_folder="/specific/netapp5/gaga/hagailevi/emp_test/true_solutions/")


