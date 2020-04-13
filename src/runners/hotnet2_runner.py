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


import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import get_network_genes


from src.runners.abstract_runner import AbstractRunner
class Hotnet2Runner(AbstractRunner):
    def __init__(self):
        super().__init__("hotnet2")

    def sif2hotnet2(self, network_file_name, script_file_name, cache_folder):

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

    def init_params(self, score_file_name, output_folder, network_file_name):

        cache_folder = os.path.join(output_folder, "cache")
        try:
            os.makedirs(cache_folder)
        except OSError:
            pass
        script_file_name=format_script(os.path.join(constants.dir_path, "src/sh", "prepare_hotnet2.sh"), ALGO_DIR=self.ALGO_DIR,
                      CACHE_DIR=cache_folder, cwd=self.ALGO_DIR)

        heat_file_name = os.path.join(cache_folder, "heatfile.txt")


        scores=pd.read_csv(score_file_name,index_col=0, sep='\t').loc[:,"qval"].apply(lambda a : -log10(a))
        scores.to_csv(heat_file_name,sep=' ',header=False)


        self.sif2hotnet2(network_file_name, script_file_name, cache_folder)
        os.remove(script_file_name)

        bg_genes=get_network_genes(network_file_name)

        return heat_file_name, bg_genes, cache_folder


    def extract_modules_and_bg(self, bg_genes, cache_folder):
        results = json.load(open(os.path.join(cache_folder, "results", "consensus", "subnetworks.json")))
        modules = [x["core"] for x in results["consensus"] if len(x["core"]) > 3 ]
        all_bg_genes = [bg_genes for x in modules]
        print("extracted {} modules".format(len(modules)))
        return modules, all_bg_genes



    def run(self, dataset_file_name, network_file_name, output_folder):
        heat_file_name, bg_genes, cache_folder = self.init_params(dataset_file_name, output_folder, network_file_name)
        script_file_name = format_script(os.path.join(constants.dir_path, "src/sh", "run_{}.sh".format(self.ALGO_NAME)),
                                         ALGO_DIR=self.ALGO_DIR,
                                         CACHE_DIR=cache_folder, OUTPUT_DIR=cache_folder,
                                         NETWORK_NAME=os.path.splitext(os.path.basename(network_file_name))[0])
        print(subprocess.Popen("bash {}".format(script_file_name), shell=True,
                               stdout=subprocess.PIPE).stdout.read())  # cwd=dir_path
        os.remove(script_file_name)
        modules, all_bg_genes = self.extract_modules_and_bg(bg_genes, cache_folder)
        print(len(modules))
        return modules, all_bg_genes


if __name__ == "__main__":
    runner=Hotnet2Runner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/dip.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}".format(runner.ALGO_NAME)))







