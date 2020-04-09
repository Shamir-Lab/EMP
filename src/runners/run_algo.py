import sys
sys.path.insert(0, '../')
import os
from multiprocessing import Process
from src import constants
from src.runners import domino_runner
from src.runners import netbox_runner
from src.runners import jactivemodules_greedy_runner
from src.runners import jactivemodules_sa_runner
# from src.runners import bionet_runner

ALGO_BY_NAMES = {"DOMINO":domino_runner.main, "netbox":netbox_runner.main, "jactivemodules_greedy":jactivemodules_greedy_runner.main, "jactivemodules_sa":jactivemodules_sa_runner.main}#, "bionet":bionet_runner.main}

def add_algo_runner(k,v):
    ALGO_BY_NAMES[k]=v

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def run_algo(dataset_name, algo, network_file_name, go_folder, output_folder, **kwargs):
    ALGO_BY_NAMES[algo](dataset_name, network_file_name, go_folder, output_folder, **kwargs)

