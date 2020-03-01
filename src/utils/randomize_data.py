import pandas as pd
import numpy as np
import shutil
import os
import src.constants as constants

import random


def get_permutation_name(prefix, dataset, algo, index):
    random_ds_name=prefix + "_random_" + dataset
    if algo is not None:
        random_ds_name += "_{}".format(algo)
    if index is not None:
        random_ds_name += "_{}".format(index)

    return random_ds_name


def permutation_output_exists(prefix, dataset, algo, index):

    return os.path.exists(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, get_permutation_name(prefix, dataset, algo, index), algo,
                     "modules_summary.tsv"))

def create_random_ds(prefix, cur_ds, index=None, algo=None):
    data_type = "score.tsv"
    cur_ds = cur_ds[len(prefix)+1:]
    random_ds_name=get_permutation_name(prefix, cur_ds, algo, index)
    root_random_dir=os.path.join(constants.DATASETS_DIR, random_ds_name)
    if os.path.exists(root_random_dir):
        shutil.rmtree(root_random_dir)
    os.makedirs(os.path.join(root_random_dir, "data"))
    os.makedirs(os.path.join(root_random_dir, "output"))
    os.makedirs(os.path.join(root_random_dir, "cache"))
    data_file_name = os.path.join(root_random_dir, "data", data_type)
    data=pd.read_csv(os.path.join(constants.DATASETS_DIR, "{}_{}".format(prefix, cur_ds), 'data', data_type), sep='\t', index_col=0)
    rd=np.random.RandomState(int(index+random.random()*10000))
    data=pd.DataFrame(data=rd.permutation(data.values), index=data.index, columns=data.columns)
    data.to_csv(data_file_name, sep='\t', index_label="id")
    return random_ds_name
