import pandas as pd
import numpy as np
import shutil
import os
import src.constants as constants

import random
from src.utils.network import get_network_genes

def get_permuted_folder_name(dataset, algo=None, index=None, ss_ratio=None):
    random_ds_name="{}_robustness".format(dataset)
    if algo is not None:
        random_ds_name += "_{}".format(algo)
    if index is not None:
        random_ds_name += "_{}".format(index)
    if ss_ratio is not None:
        random_ds_name += "_{}".format(ss_ratio)

    return random_ds_name


def permutation_solution_exists(dataset, algo, index, ss_ratio, output_folder):

    return os.path.exists(
        os.path.join(output_folder, "sol_{}_{}".format(algo, get_permuted_folder_name(dataset, index, ss_ratio)), "report",
                     "modules_summary.tsv"))

def permutation_dataset_exists(dataset, index, ss_ratio, output_folder):

    return os.path.exists(
        os.path.join(output_folder, get_permuted_folder_name(dataset, index, ss_ratio), "data", "scores.tsv"))

def create_random_ds(output_folder,dataset_file, network_file_name, ss_ratio, index=None, algo=None):
    dataset=os.path.splitext(os.path.split(dataset_file)[1])[0]
    random_ds_name=get_permuted_folder_name(dataset, algo, index, ss_ratio)
    df_scores=pd.read_csv(os.path.join(dataset_file), sep='\t', index_col=0)
    assay_genes = set(df_scores.index.values)
    # network_genes = get_network_genes(network_file_name)
    # overlapped_genes = list(network_genes.intersection(assay_genes))
    overlapped_genes=assay_genes
    rd=np.random.RandomState(int(index+random.random()*10000))
    # df_permuted_scores=pd.DataFrame(data=rd.permutation(df_scores.values), index=df_scores.index, columns=df_scores.columns)
    genes_to_omit = np.random.choice(np.array(list(overlapped_genes)), int(len(overlapped_genes) * ss_ratio), replace=False)
    df_permuted_scores = df_scores.loc[~df_scores.index.isin(genes_to_omit), :]
    premuted_scores_folder=os.path.join(output_folder, random_ds_name)
    if os.path.exists(premuted_scores_folder):
        shutil.rmtree(premuted_scores_folder)
    os.makedirs(os.path.join(premuted_scores_folder,"data"))
    os.makedirs(os.path.join(premuted_scores_folder,"output"))
    permuted_scored_file = os.path.join(premuted_scores_folder, "data",constants.SCORES_FILE_NANE)
    df_permuted_scores.to_csv(permuted_scored_file, sep='\t', index_label="id")
    return premuted_scores_folder
