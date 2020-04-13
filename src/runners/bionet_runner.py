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

from src.utils.network import get_network_genes
from src.utils.network import remove_subgraph_by_nodes
import numpy as np
import src.constants as constants

from src.runners.abstract_runner import AbstractRunner
class BionetRunner(AbstractRunner):
    def __init__(self):
        super().__init__("bionet")



def run_bionet(self, deg_file_name, network_file_name, fdr=0.05):
    from src.utils.r_runner import run_rscript
    script = open(os.path.join(self.ALGO_DIR, "{}.r".format(self.ALGO_NAME))).read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name, fdr=fdr)


def init_params(self, network_file_name, omitted_genes = [], ts=str(time.time())):
    new_network_file_name=remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)
    bg_genes=get_network_genes(network_file_name)

    return new_network_file_name, bg_genes

def get_module(self, network_file_name, score_file_name, omitted_genes, ts=str(time.time()),fdr=0.05):
    network_file_name, bg_genes = init_params(network_file_name=network_file_name, omitted_genes=omitted_genes)
    results = run_bionet(score_file_name, network_file_name,fdr)
    module_genes = np.array(results["module_genes"])

    
    sys.stdout.write("module gene size: {}. ratio: {}\n".format(module_genes.shape[0], module_genes.shape[0]/float(len(bg_genes))))
    return list(module_genes), bg_genes


def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):

    fdr=0.05

    omitted_genes = []
    modules = []
    all_bg_genes = []
    small_modules=0
    for x in range(50):
        try:
           module_genes, bg_genes = get_module(network_file_name, dataset_file_name, omitted_genes, str(x), fdr=fdr)
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
    return modules, all_bg_genes



if __name__ == "__main__":
    runner=BionetRunner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/dip.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}".format(runner.ALGO_NAME)))









