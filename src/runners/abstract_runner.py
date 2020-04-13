#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../..')
import os
import src.constants as constants
from src.utils.network import build_all_reports

class AbstractRunner(object):

    def __init__(self, algo_name):
        self.ALGO_NAME=algo_name
        self.ALGO_DIR = os.path.join(constants.ALGO_DIR, self.ALGO_NAME)


    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        """implement. return modules, all_bg_genes"""

    def build_all_reports(self, algo_name, modules, all_bg_genes, go_folder, output_folder):
        build_all_reports(algo_name, modules, all_bg_genes, go_folder, output_folder)

    def main(self, dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):
        modules, all_bg_genes = self.run(os.path.abspath(dataset_file_name), os.path.abspath(network_file_name), os.path.abspath(output_folder), **kwargs)
        self.build_all_reports(self.ALGO_NAME, modules, all_bg_genes, os.path.abspath(go_folder), os.path.abspath(os.path.join(output_folder, "report")))







