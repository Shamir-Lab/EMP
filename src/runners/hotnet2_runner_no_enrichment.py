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
from numpy import log10
import pandas as pd
import subprocess
import simplejson as json
import shutil

import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import get_network_genes

ALGO_NAME = "hotnet2"
ALGO_DIR = os.path.join(constants.ALGO_DIR, ALGO_NAME)

def sif2hotnet2(network_file_name, script_file_name, cache_folder):

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))

    lns = ["{} {}".format(i+1, cur) for i, cur in enumerate(vertices)]
    open(os.path.join(cache_folder, "hotnet2_vertices.txt"), "w+").write("\n".join(lns))

    inds = map(lambda x: (vertices.index(x[0])+1, vertices.index(x[1])+1), zip(src,dst))


    lns = ["{} {} {}".format(cur[0], cur[1], 1) for cur in inds]
    open(os.path.join(cache_folder, "hotnet2_edges.txt"), "w+").write("\n".join(lns))

    print(subprocess.Popen("bash {}".format(script_file_name) , shell=True, stdout=subprocess.PIPE).stdout.read())

def init_params(score_file_name, output_folder, network_file_name):

    cache_folder = os.path.join(output_folder, "cache")
    try:
        os.makedirs(cache_folder)
    except OSError:
        pass
    script_file_name=format_script(os.path.join(constants.dir_path, "src/sh", "prepare_hotnet2.sh"), ALGO_DIR=ALGO_DIR,
                  CACHE_DIR=cache_folder, cwd=ALGO_DIR)

    heat_file_name = os.path.join(cache_folder, "heatfile.txt")


    scores=pd.read_csv(score_file_name,index_col=0, sep='\t').loc[:,"qval"].apply(lambda a : -log10(a) if a != 0 else 315)
    scores.to_csv(heat_file_name,sep=' ',header=False)


    # import src.utils.infra as infra
    # deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
    # h_rows, h_cols, deg_data = infra.separate_headers(deg)
    # ind = np.where(h_cols == "qval")[0][0]

    # lns = []
    # for i, cur in enumerate(deg_data):
    #     lns.append(" ".join([str(h_rows[i]), str(-log10(cur[ind]))]))
    # open(heat_file_name, "w+").write("\n".join(lns))


    sif2hotnet2(network_file_name, script_file_name, cache_folder)
    os.remove(script_file_name)

    bg_genes=get_network_genes(network_file_name)

    return heat_file_name, bg_genes, cache_folder


def extract_modules_and_bg(bg_genes, cache_folder):
    results = json.load(open(os.path.join(cache_folder, "results", "consensus", "subnetworks.json")))
    modules = [x["core"] for x in results["consensus"] if len(x["core"]) > 3 ]
    all_bg_genes = [bg_genes for x in modules]
    print("extracted {} modules".format(len(modules)))
    shutil.rmtree(cache_folder)
    os.makedirs(cache_folder)
    for i, genes in enumerate(modules):
        modules_genes_file_name=os.path.join(cache_folder, "hotnet2_module_{}_genes.txt".format(i))
        open(modules_genes_file_name, "w+").write("\n".join(genes))
    return modules, all_bg_genes


def main(dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):

    heat_file_name, bg_genes, cache_folder = init_params(dataset_file_name, output_folder, network_file_name)

    script_file_name=format_script(os.path.join(constants.dir_path, "src/sh", "run_{}.sh".format(ALGO_NAME)), ALGO_DIR=ALGO_DIR,
                  CACHE_DIR=cache_folder, OUTPUT_DIR=cache_folder, NETWORK_NAME=os.path.splitext(os.path.basename(network_file_name))[0])
    print(subprocess.Popen("bash {}".format(script_file_name), shell=True,
                           stdout=subprocess.PIPE).stdout.read())  # cwd=dir_path
    # os.remove(script_file_name)
    modules, all_bg_genes = extract_modules_and_bg(bg_genes, cache_folder)
    print(len(modules))
    # build_all_reports(ALGO_NAME, modules, all_bg_genes, go_folder, os.path.join(output_folder, "report"))


if __name__ == "__main__":

    main(dataset_file_name="/media/hag007/Data/emp_test/original_datasets/brca.tsv", network_file_name="/media/hag007/Data/emp_test/networks/dip.sif", go_folder="/media/hag007/Data/emp_test/go", output_folder="/media/hag007/Data/emp_test/true_solutions/brca_hotnet2")







