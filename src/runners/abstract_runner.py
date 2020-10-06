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
import time
from src.utils.network import build_all_reports
import signal


class AbstractRunner(object):

    def __init__(self, algo_name):
        self.ALGO_NAME=algo_name
        self.ALGO_DIR = os.path.join(constants.ALGO_DIR, self.ALGO_NAME)

    def handler(self, signum, frame):
        print("Forever is over!")
        try:
           os.makedirs(os.path.split(self.full_path)[0])
        except Exception as e:
           print(e)
     
        open(self.full_path, 'w+').write("9999")
 
        raise Exception("end of time")

    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        """implement. return modules, all_bg_genes"""

    def build_all_reports(self, algo_name, modules, all_bg_genes, go_folder, output_folder):
        build_all_reports(algo_name, modules, all_bg_genes, go_folder, output_folder)

    def main_timer(self, dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):
        print(f'about to run algo with the following parameters: {os.path.abspath(dataset_file_name)}, {os.path.abspath(network_file_name)}, {os.path.abspath(output_folder)}, {kwargs}')

        basic_name=output_folder.split('/')[-1]
        print("basic_name: {}".format(basic_name))
        if "random" in basic_name:

            self.full_path=os.path.join(constants.TIMER_DIR, os.path.splitext(os.path.basename(network_file_name))[0], basic_name.split("_")[-3], "_".join(basic_name.split("_")[1:-3]+ basic_name.split("_")[-1:]))+".txt" 
        else:
            self.full_path=os.path.join(constants.TIMER_DIR, os.path.splitext(os.path.basename(network_file_name))[0], basic_name.split("_")[0], "_".join(basic_name.split("_")[1:]))+".txt"

        signal.signal(signal.SIGALRM, self.handler)
        signal.alarm(int(60*60*100))
        tic = time.perf_counter() 
        try:
            modules, all_bg_genes = self.run(os.path.abspath(dataset_file_name), os.path.abspath(network_file_name), os.path.abspath(output_folder), **kwargs)
            toc = time.perf_counter()
        except Exception as e:
            if str(e)=="end of time":
                return
            else:
                raise e
               
        print(f"algo execution time: {toc-tic:0.2f}")
 
        try:
           os.makedirs(os.path.split(self.full_path)[0])
        except Exception as e:
           print(e)
     
        open(self.full_path, 'w+').write(f"{toc-tic:0.2f}")
  
        print(f"calculate GO enrichment for {len(modules)} modules")
        self.build_all_reports(self.ALGO_NAME, modules, all_bg_genes, os.path.abspath(go_folder), os.path.abspath(os.path.join(output_folder, "report")))



    def main(self, dataset_file_name, network_file_name, go_folder, output_folder, **kwargs):
        print(f'about to run algo with the following parameters: {os.path.abspath(dataset_file_name)}, {os.path.abspath(network_file_name)}, {os.path.abspath(output_folder)}, {kwargs}')

        modules, all_bg_genes = self.run(os.path.abspath(dataset_file_name), os.path.abspath(network_file_name), os.path.abspath(output_folder), **kwargs)
