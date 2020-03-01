import sys
sys.path.insert(0, '../')

import argparse
import src.constants as constants
import os
import time
import shutil
import pandas as pd
from src.utils.network import get_network_genes

def get_score(dataset_name):
    score_file_name = os.path.join(constants.DATASETS_DIR,dataset_name,"data", "score.tsv".format(score_method))
    df = pd.read_csv(score_file_name, sep="\t")
    if not "pval" in df or not  "qval" in df:
        raise Exception("expected cols: pval, qval or score. got {}".format(df.columns))

    return score_file_name

def init_common_params(NETWORK_NAME, dataset_name):

    network_file_name = os.path.join(constants.NETWORKS_DIR, NETWORK_NAME)
    score_file_name = get_score(dataset_name)
    network_genes = get_network_genes()
    return network_file_name, score_file_name, network_genes

