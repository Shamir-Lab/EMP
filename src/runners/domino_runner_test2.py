import sys
sys.path.insert(0, '../')

import os
import pandas as pd

from src import constants
from src.implementations.domino import main as domino_main
from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.network import get_network_genes

from src.runners.abstract_runner import AbstractRunner
class DominoRunner(AbstractRunner):
    def __init__(self):
        super().__init__("DOMINO_test2")


    def extract_modules_and_bg(self, bg_genes, dest_algo_dir):
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

    def init_params(self, dataset_file_name, network_file_name, output_folder):


        df_scores=pd.read_csv(dataset_file_name, sep='\t', index_col=0)
        sig_genes=df_scores['qval'][df_scores['qval']<0.05].index
        active_genes_file=os.path.join(output_folder, "active_genes_file.txt")
        open(active_genes_file, "w+").write("\n".join([x for x in sig_genes if len(ensembl2entrez_convertor([x]))>0 ]))
        bg_genes=get_network_genes(network_file_name)
        return active_genes_file, bg_genes


    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        print("run domino runner...")
        slices_file = kwargs['slices_file']
        constants.N_OF_THREADS=1
        if 'n_of_threads' in kwargs:
            constants.N_OF_THREADS=kwargs['n_of_threads']
        constants.USE_CACHE=False
        if 'use_cache' in kwargs:
            constants.USE_CACHE=kwargs['use_cache']=='true'
        slice_threshold = 0.3
        if 'slice_threshold' in kwargs:
            slice_threshold = kwargs['slice_threshold']
        module_threshold = 0.05
        if 'module_threshold' in kwargs:
            module_threshold = kwargs['module_threshold']
        active_genes_file, bg_genes = self.init_params(dataset_file_name, network_file_name, output_folder)
        print(f'domino_parameters: active_genes_file={active_genes_file}, network_file={network_file_name},slices_file={slices_file}, slice_threshold={slice_threshold},module_threshold={module_threshold}')
        modules = domino_main(active_genes_file=active_genes_file, network_file=network_file_name,
                              slices_file=slices_file, slice_threshold=slice_threshold,
                              module_threshold=module_threshold)
        modules = list(filter(lambda x: len(x) > 3, modules))
        all_bg_genes = [bg_genes for x in modules]
        return modules, all_bg_genes


if __name__ == "__main__":
    constants.N_OF_THREADS=1
    constants.USE_CACHE=False
    runner=DominoRunner()
    runner.main(dataset_file_name=os.path.join(constants.config_json["base_dir"],"original_datasets/brca.tsv"), network_file_name=os.path.join(constants.config_json["base_dir"],"networks/string_minimal_no_prefix_500_g.sif"), go_folder=os.path.join(constants.config_json["base_dir"],"go"), output_folder=os.path.join(constants.config_json["base_dir"],"true_solutions/tnfa_{}".format(runner.ALGO_NAME)), slices_file=os.path.join(constants.config_json["base_dir"],"networks/string_modules_agg.txt"), module_threshold=0.05, slice_threshold=0.3)




