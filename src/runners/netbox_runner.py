#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys
sys.path.insert(0, '../')


import os
import subprocess
import pandas as pd
import shutil

from src import constants

from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.scripts import format_script
from src.utils.network import get_network_genes

import random

from src.runners.abstract_runner import AbstractRunner

class NetboxRunner(AbstractRunner):

    def __init__(self):
        super().__init__("netbox")

    def init_params(self, dataset_file_name, network_file_name, output_folder, dest_algo_dir):

        df_scores=pd.read_csv(dataset_file_name, sep='\t', index_col=0)
        sig_genes=df_scores['qval'][df_scores['qval']<0.05].index
        active_genes_file=os.path.join(output_folder, "active_genes_file.txt")
        open(active_genes_file, "w+").write("\n".join([x for x in sig_genes if len(ensembl2entrez_convertor([x]))>0 ]))
        bg_genes=get_network_genes(network_file_name)
        conf_file = "conf.props"
        conf_file_name=format_script(os.path.join(dest_algo_dir, conf_file), pval_threshold=0.05, sp_threshold=2, gene_file=active_genes_file)
        return active_genes_file, bg_genes, conf_file_name

    def extract_modules_and_bg(self, dest_algo_dir):
        results = open(os.path.join(dest_algo_dir, "modules.txt")).readlines()
        modules = [[] for x in range(max([0]+[int(x.strip().split(" =")[1]) for x in results[1:]]) + 1)]
        for x in results[1:]:
            if int(x.strip().split(" =")[1]) != -1:
                modules[int(x.strip().split(" =")[1])].append(x.strip().split(" =")[0])
            else:
                modules.append([x.strip().split(" =")[0]])
        modules = list(filter(lambda x: len(x) > 3, modules))
        print("extracted {} modules".format(len(modules)))
        return modules


    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        script_name = "run_{}.sh".format(self.ALGO_NAME)
        ######### client (single-threaded)
        # dest_algo_dir="{}".format(ALGO_DIR,random.random())
        # conf_file_name, bg_genes, conf_file_name=init_params(dataset_file_name, network_file_name, output_folder, dest_algo_dir)
        ######### server (multi-threaded)
        dest_algo_dir = "{}_{}".format(self.ALGO_DIR, random.random())
        shutil.copytree(self.ALGO_DIR, dest_algo_dir)
        conf_file_name, bg_genes, conf_file_name = self.init_params(dataset_file_name, network_file_name, output_folder,
                                                               dest_algo_dir)
        #########
        try:
            os.remove(os.path.join(dest_algo_dir, "modules.txt"))
        except:
            pass
        script_file_name = format_script(os.path.join(constants.dir_path, "src/sh", script_name),
                                         NETBOX_HOME=dest_algo_dir,
                                         BASE_FOLDER=os.path.abspath("../../data/emp_test"),
                                         CONFIG_FILE_NAME=conf_file_name)
        print(subprocess.Popen("bash {}".format(script_file_name), shell=True,
                               stdout=subprocess.PIPE, cwd=dest_algo_dir).stdout.read())
        os.remove(script_file_name)
        modules = self.extract_modules_and_bg(dest_algo_dir)
        all_bg_genes = [bg_genes for x in modules]
        ######### server (on parallel)
        shutil.rmtree(dest_algo_dir)

        return modules, all_bg_genes




if __name__ == "__main__":
    runner=NetboxRunner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/dip.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}".format(runner.ALGO_NAME)))






