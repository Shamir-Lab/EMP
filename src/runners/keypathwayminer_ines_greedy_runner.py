#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../..')

import os
import subprocess

import shutil
import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import build_all_reports
from src.utils.network import remove_subgraph_by_nodes

from src.utils.network import get_network_genes

import pandas as pd


import random
from src.runners.abstract_runner import AbstractRunner
class KPMRunner(AbstractRunner):
    def __init__(self):
        super().__init__("keypathwayminer")
        self.strategy = "INES"
        self.algorithm = "GREEDY"

    def init_params(self, score_file_name, omitted_genes, network_file_name, ts, dest_algo_dir):
        if os.path.exists(os.path.join(dest_algo_dir, "results")):
            shutil.rmtree(os.path.join(dest_algo_dir, "results"))

        binary_score_file_name = os.path.join(dest_algo_dir, "binary_score.txt")

        df=pd.read_csv(score_file_name, sep='\t', index_col=0).loc[:,"qval"].sort_values()
        df=(df<=0.05).astype(int)
        df.to_csv(binary_score_file_name, sep='\t')

        new_network_file_name = remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)
        bg_genes=get_network_genes(network_file_name)

        return binary_score_file_name, new_network_file_name, bg_genes


    def format_scripts(self, score_file_name, network_name="dip", STRATEGY="INES", algorithm="GREEDY", algo_dir=None, dataset_name=None):
        script_file_name=format_script(os.path.join(constants.dir_path,"src","sh", "run_{}.sh".format(self.ALGO_NAME)), STRATEGY=STRATEGY, ALGORITHM=algorithm)
        format_script(os.path.join(algo_dir, "kpm.properties"), uid=False, network_name=network_name, algorithm=algorithm)
        format_script(os.path.join(algo_dir, "datasets_file.txt"), uid=False, score_file_name=score_file_name)

        return script_file_name


    def extract_module_genes(self, bg_genes, dest_algo_dir):
        i = 1
        modules = []
        while os.path.exists(os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))) and i<2:
            results = open(
                os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))).readlines()
            results = list(map(lambda x: x.strip(), results))
            modules.append(results)
            i += 1

        return modules, [bg_genes for x in modules]



    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        omitted_genes = []
        modules = []
        all_bg_genes = []
        dest_algo_dir = "{}_{}".format(self.ALGO_DIR, random.random())
        shutil.copytree(self.ALGO_DIR, dest_algo_dir)
        empty_counter = 0
        for cur_i_module in range(50):
            binary_score_file_name, cur_network_file_name, bg_genes = self.init_params(dataset_file_name, omitted_genes,
                                                                                       network_file_name,
                                                                                       str(random.random()),
                                                                                       dest_algo_dir)

            script_file_name = self.format_scripts(score_file_name=binary_score_file_name,
                                                   network_name=cur_network_file_name,
                                                   STRATEGY=self.strategy, algorithm=self.algorithm, algo_dir=dest_algo_dir)
            print(subprocess.Popen("bash {}".format(script_file_name), shell=True,
                                   stdout=subprocess.PIPE, cwd=dest_algo_dir).stdout.read())
            module, all_bg_gene = self.extract_module_genes(bg_genes, dest_algo_dir)
            print(module)
            if len(module[0]) > 3:
                empty_counter = 0
                modules.append(module[0])
                all_bg_genes.append(all_bg_gene[0])
            else:
                empty_counter += 1
            omitted_genes += list(module[0])
            os.remove(script_file_name)

            if empty_counter > 3:
                print("got more that 3 small modules in row. continue...")
                break
        shutil.rmtree(dest_algo_dir)
        return modules, all_bg_genes


if __name__ == "__main__":
    runner=KPMRunner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/dip.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}_{}_{}".format(runner.ALGO_NAME, runner.strategy, runner.algorithm)))

