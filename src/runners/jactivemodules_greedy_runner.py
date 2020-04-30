"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../')


import os
import subprocess

import src.constants as constants

from src.utils.scripts import format_script
from src.utils.network import get_network_genes

from src.runners.abstract_runner import AbstractRunner
class jAMGreedyRunner(AbstractRunner):
    def __init__(self):
        super().__init__("jactivemodules")
        self.search_method="greedy"

    def init_params(self, network_file_name, output_folder):
        bg_genes=get_network_genes(network_file_name)
        results_file_name = "{}/{}_{}_results.txt".format(output_folder, self.ALGO_NAME, self.search_method)
        return results_file_name, bg_genes


    def extract_modules_and_bg(self, bg_genes, results_file_name, modules_genes_file_name, output_folder):
        results = open(results_file_name.format(output_folder, self.ALGO_NAME)).readlines()
        modules = [x.split()[:-1] for x in results]
        modules = [cur for cur in modules if len(cur) > 3]
        all_bg_genes = [bg_genes for x in modules]
        module_genes = [y for x in modules for y in x]
        open(modules_genes_file_name, "w+").write("\n".join(module_genes))
        print("extracted {} modules".format(len(modules)))
        return all_bg_genes, modules



    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        print("start running jactivemodules_greedy")
        results_file_name, bg_genes = self.init_params(network_file_name, output_folder)
        script_file_name = format_script(os.path.join(constants.dir_path, "src/sh", "run_{}.sh".format(self.ALGO_NAME)),
                                         BASE_FOLDER="/specific/netapp5/gaga/hagailevi/emp_test",
                                         ALGO_DIR=self.ALGO_DIR, NETWORK_FILE_NAME=network_file_name,
                                         SCORE_FILE_NAME=dataset_file_name,
                                         IS_GREEDY=str(self.search_method == "greedy"), OUTPUT_FILE=results_file_name,
                                         NUM_OF_MODULES=50, OVERLAP_THRESHOLD=0)
        print("running :{}".format(script_file_name))
        subprocess.Popen("bash {}".format(script_file_name), shell=True,
                         stdout=subprocess.PIPE, cwd=self.ALGO_DIR).stdout.read()
        os.remove(script_file_name)
        modules_genes_file_name = os.path.join(output_folder,
                                               "{}_{}_module_genes.txt".format(self.ALGO_NAME, self.search_method))

        all_bg_genes, modules = self.extract_modules_and_bg(bg_genes, results_file_name, modules_genes_file_name,
                                                            output_folder)
        print([len(m) for m in modules])
        return modules, all_bg_genes


if __name__ == "__main__":
    from src.utils.go import init_state
    init_state(os.path.join(constants.config_json["base_dir"],"go"))
    runner=jAMGreedyRunner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/tnfa.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/dip.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}_{}".format(runner.ALGO_NAME, runner.search_method)))


