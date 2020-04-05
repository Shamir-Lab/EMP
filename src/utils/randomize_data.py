import pandas as pd
import numpy as np
import shutil
import os
import src.constants as constants

import random


def get_permuted_folder_name(dataset, algo=None, index=None):
    random_ds_name="{}_random".format(dataset)
    if algo is not None:
        random_ds_name += "_{}".format(algo)
    if index is not None:
        random_ds_name += "_{}".format(index)

    return random_ds_name


def permutation_solution_exists(dataset, algo, index, output_folder):

    return os.path.exists(
        os.path.join(output_folder, "sol_{}_{}".format(algo, get_permuted_folder_name(dataset, index)), "report",
                     "modules_summary.tsv"))

def permutation_dataset_exists(dataset, index, output_folder):
    
    return os.path.exists(
        os.path.join(output_folder, get_permuted_folder_name(dataset, index), "data", "scores.tsv"))

def create_random_ds(output_folder,dataset_file, index=None, algo=None):
    dataset=os.path.splitext(os.path.split(dataset_file)[1])[0]
    random_ds_name=get_permuted_folder_name(dataset, algo, index)
    df_scores=pd.read_csv(os.path.join(dataset_file), sep='\t', index_col=0)
    rd=np.random.RandomState(int(index+random.random()*10000))
    df_permuted_scores=pd.DataFrame(data=rd.permutation(df_scores.values), index=df_scores.index, columns=df_scores.columns)

    premuted_scores_folder=os.path.join(output_folder, random_ds_name)
    if os.path.exists(premuted_scores_folder):
        shutil.rmtree(premuted_scores_folder)
    os.makedirs(os.path.join(premuted_scores_folder,"data"))
    os.makedirs(os.path.join(premuted_scores_folder,"output"))
    permuted_scored_file = os.path.join(premuted_scores_folder, "data",constants.SCORES_FILE_NANE)
    df_permuted_scores.to_csv(permuted_scored_file, sep='\t', index_label="id")
    return premuted_scores_folder
