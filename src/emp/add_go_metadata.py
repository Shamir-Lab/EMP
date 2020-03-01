import sys
sys.path.insert(0, '../')

import src.utils.add_GO_terms_metadata_agg
import pandas as pd
import os
from src import constants
import argparse

def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset', dest='dataset', default="TNFa_2")
    parser.add_argument('--algo', dest='algo', default="my_netbox_td")
    parser.add_argument('--omic_type', dest='omic_type', default="GE")
    parser.add_argument('--output_folder', dest='output_folder', default=None)
    parser.add_argument('--n_permutations', dest='n_permutations', default=5000)

    args = parser.parse_args()

    dataset=args.dataset.split(",")
    algo=args.algo.split(",")
    omic_type = args.omic_type
    output_folder = args.output_folder
    if output_folder is None:
        output_folder=os.path.join(dataset,"output")
    n_permutations=int(args.n_permutations)

    print("{}_{}".format(dataset, algo))
    csv_file_name=os.path.join(output_folder, "MAX/emp_diff_modules_{dataset}_{algo}.tsv")
    df_all=src.utils.add_GO_terms_metadata_agg.add_md_to_terms(dataset, algo, n_permutations, csv_file_name=csv_file_name, omic_type=omic_type)
    df_all.loc[df_all["hg_pval"].values > 0, :][
        ["GO name", "hg_pval", "hg_pval_max", "emp_pval", "hg_rank", "emp_rank", "n_genes", "depth",
         "passed_fdr"]].to_csv(csv_file_name[:-4] + "_md.tsv", sep='\t')


if __name__=="__main__":
    main()