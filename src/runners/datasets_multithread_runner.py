import sys
sys.path.insert(0, '../')
import os
from multiprocessing import Process
from src import constants
from src.runners import domino_runner

ALGO_BY_NAMES = {"domino":domino_runner.main}

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def run_dataset(dataset_name, expected_genes=None, disease_name=None, score_method=constants.DEG_EDGER, algos=None, output_folder=None, network_file_name="dip.sif"):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    prcs = []
    for cur in algos:
        prcs.append(Process(target=ALGO_BY_NAMES[cur], args=[dataset_name, disease_name, expected_genes, score_method, output_folder, network_file_name]))
        prcs[-1].start()

    for cur in prcs:
        cur.join()


if __name__ == "__main__":

    datasets=["EN_CML", "EN_LICH", "EN_BRCA", "EN_KIRC"] # ["EN_LUNG", "EN_PRAD", "EN_PAAD"]
    algos = ["hotnet2", "bionet", "jactivemodules_greedy", "jactivemodules_sa", "netbox", "keypathwayminer_INES_GREEDY"]

    for cur_ds in datasets:
        print("current folder : {}".format(os.path.basename(cur_ds)))
        score_method = constants.PREDEFINED_SCORE

        run_dataset(cur_ds, score_method=score_method,
                    algos=algos, network_file_name="dip.sif") #
