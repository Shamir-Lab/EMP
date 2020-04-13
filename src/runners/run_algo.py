import sys
sys.path.insert(0, '../')
import os
from src import constants
# from src.runners import domino_runner
from src.runners.netbox_runner import NetboxRunner
from src.runners.jactivemodules_greedy_runner import jAMGreedyRunner
from src.runners.jactivemodules_sa_runner import jAMSARunner
from src.runners.bionet_runner import BionetRunner
from src.runners.keypathwayminer_ines_greedy_runner import KPMRunner

ALGO_BY_NAMES = {"netbox":NetboxRunner(), "jactivemodules_greedy": jAMGreedyRunner(), "jactivemodules_sa": jAMSARunner(), "bionet": BionetRunner(), 'keypathwayminer_INES_GREEDY': KPMRunner()}


def add_algo_runner(k,v):
    ALGO_BY_NAMES[k]=v

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def run_algo(dataset_name, algo, network_file_name, go_folder, output_folder, **kwargs):
    ALGO_BY_NAMES[algo].main(dataset_name, network_file_name, go_folder, output_folder, **kwargs)


