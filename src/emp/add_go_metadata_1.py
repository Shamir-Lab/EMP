import sys
sys.path.insert(0, '../..')


import pandas as pd
import os
from src import constants
import argparse

def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data1/emp_test/datasets/scz.tsv")
    parser.add_argument('--algo', dest='algo', default="DOMINO")
    parser.add_argument('--go_folder', dest='go_folder', default="/media/hag007/Data1/emp_test/go")
    parser.add_argument('--report_folder', dest='report_folder', default="/media/hag007/Data1/emp_test/report")
    parser.add_argument('--n_permutations', dest='n_permutations', default=10 )

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    report_folder=args.report_folder
    go_folder=args.go_folder
    n_permutations=int(args.n_permutations)

    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
    constants.GO_DIR=go_folder
    import src.utils.add_GO_terms_metadata
    print("{}_{}".format(dataset_name, algo))

    csv_file_name=os.path.join(report_folder, "emp_diff_modules_{dataset}_{algo}.tsv")

    df_all=src.utils.add_GO_terms_metadata.add_md_to_terms(dataset_name, algo, n_permutations, csv_file_name=csv_file_name)
    df_all.loc[:, :][
        ["GO name", "hg_pval", "hg_pval_max", "hg_rank", "n_genes", "depth"]].to_csv(csv_file_name.format(dataset=dataset_name, algo=algo)[:-4] + "_md.tsv", sep='\t') # df_all["hg_pval_max"].values > 0 # "passed_fdr", "emp_pval", "emp_rank"


if __name__=="__main__":
    main()