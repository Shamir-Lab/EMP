"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../')


import os
import numpy as np
import pandas as pd
import subprocess
import json


import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import build_all_reports, get_network_genes
import src.utils.server as server
from src.utils.network import output_modules

ALGO_NAME = "jactivemodules"
ALGO_DIR = os.path.join(constants.ALGO_DIR, ALGO_NAME)


def init_params(search_method, network_file_name, output_folder):
    bg_genes=get_network_genes(network_file_name)
    results_file_name = "{}/{}_{}_results.txt".format(output_folder, ALGO_NAME, search_method)
    return results_file_name, bg_genes


def extract_modules_and_bg(bg_genes, results_file_name, modules_genes_file_name, output_folder):
    results = open(results_file_name.format(output_folder, ALGO_NAME)).readlines()
    modules = [x.split()[:-1] for x in results]
    modules = [cur for cur in modules if len(cur) > 3]
    all_bg_genes = [bg_genes for x in modules]
    module_genes = [y for x in modules for y in x]
    open(modules_genes_file_name, "w+").write("\n".join(module_genes))
    print("extracted {} modules".format(len(modules)))
    return all_bg_genes, modules


def main(dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):
    print("start running jactivemodules_greedy")
    search_method = "greedy"
    results_file_name, bg_genes = init_params(search_method, network_file_name, output_folder)


    script_file_name=format_script(os.path.join(constants.dir_path, "src/sh", "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER="/specific/netapp5/gaga/hagailevi/emp_test",
                  ALGO_DIR=ALGO_DIR, NETWORK_FILE_NAME=network_file_name, SCORE_FILE_NAME=dataset_file_name,
                  IS_GREEDY=str(search_method == "greedy"), OUTPUT_FILE=results_file_name, NUM_OF_MODULES=50, OVERLAP_THRESHOLD=0)
    print("running :{}".format(script_file_name))
    subprocess.Popen("bash {}".format(script_file_name), shell=True,
                     stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    os.remove(script_file_name)
    modules_genes_file_name = os.path.join(output_folder, "{}_{}_module_genes.txt".format(ALGO_NAME, search_method))
    all_bg_genes, modules = extract_modules_and_bg(bg_genes, results_file_name, modules_genes_file_name, output_folder)

    build_all_reports(ALGO_NAME, modules, all_bg_genes, go_folder, os.path.join(output_folder, "report"))


if __name__ == "__main__":

    main(dataset_file_name=os.path.join(constants.config_json,"emp_test/original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json,"emp_test/networks/dip.sif"), go_folder=os.path.join(constants.config_json,"emp_test/networks/dip.sif"), output_folder=os.path.join(constants.config_json,"hagailevi/emp_test/true_solutions/tnfa_{}".format(ALGO_NAME)))


