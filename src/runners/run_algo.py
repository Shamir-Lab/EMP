import sys
sys.path.insert(0, '../')
import os
from src import constants
from src.runners.domino_runner import DominoRunner
from src.runners.domino_runner_string import DominoRunner as DominoRunnerString
from src.runners.domino_runner_test import DominoRunner as DominoRunnerTest
from src.runners.domino_runner_huri import DominoRunner as DominoRunnerHuri
from src.runners.domino_runner_test2 import DominoRunner as DominoRunnerTest2
from src.runners.domino_runner_test3 import DominoRunner as DominoRunnerTest3
from src.runners.domino_runner_test4 import DominoRunner as DominoRunnerTest4
from src.runners.domino_runner_test5 import DominoRunner as DominoRunnerTest5
from src.runners.domino_runner_test6 import DominoRunner as DominoRunnerTest6
from src.runners.netbox_runner import NetboxRunner
from src.runners.jactivemodules_greedy_runner import jAMGreedyRunner
from src.runners.jactivemodules_sa_runner import jAMSARunner
from src.runners.bionet_runner import BionetRunner
from src.runners.keypathwayminer_ines_greedy_runner import KPMRunner
from src.runners.hotnet2_runner import Hotnet2Runner

ALGO_BY_NAMES = {"netbox":NetboxRunner(), "netbox2_string":NetboxRunner(), "netbox2":NetboxRunner(), "netbox3":NetboxRunner(), "jactivemodules_greedy": jAMGreedyRunner(), "jactivemodules_greedy_string": jAMGreedyRunner(), "jactivemodules_sa": jAMSARunner(), "bionet": BionetRunner(), "bionet_string": BionetRunner(), 'keypathwayminer_INES_GREEDY': KPMRunner(), 'hotnet2': Hotnet2Runner(), 'DOMINO2': DominoRunner(), 'DOMINO3': DominoRunnerString(), 'DOMINO4' : DominoRunnerHuri(), 'DOMINO_test' : DominoRunnerTest(), 'DOMINO_test2' : DominoRunnerTest2(), 'DOMINO_test3' : DominoRunnerTest3(), 'DOMINO_test4' : DominoRunnerTest4(), 'DOMINO_test5' : DominoRunnerTest5(), 'DOMINO_test6' : DominoRunnerTest6()}


def add_algo_runner(k,v):
    ALGO_BY_NAMES[k]=v

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def run_algo(dataset_name, algo, network_file_name, go_folder, output_folder, **kwargs):
    ALGO_BY_NAMES[algo].main(dataset_name, network_file_name, go_folder, output_folder, **kwargs)


