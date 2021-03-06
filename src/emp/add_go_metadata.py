import sys
sys.path.insert(0, '../..')
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from functools import reduce
from src.utils.go import get_all_genes_for_term
import pandas as pd
import os
import src.constants as constants
import argparse
from src.utils.go import init_state

def add_md_to_terms(dataset="SOC", algo="jactivemodules_sa", n_permutations=300, csv_file_name=None, filtered_go_ids_file="", hg_th=0.001):

    csv_file_name = csv_file_name.format(dataset=dataset, algo=algo)
    try:
        df = pd.read_csv(csv_file_name, sep='\t', index_col=0)
    except:
        print("no bg file {} for {}, {}".format(csv_file_name, dataset,algo))
        return None
    df = df.dropna()

    n_genes = [len(get_all_genes_for_term(cur_go_id, cur_go_id, cur_go_id == cur_go_id)) for i, cur_go_id in
               enumerate(df.index.values)]
    df["n_genes"] = pd.Series(n_genes, index=df.index)
    df["hg_pval_max"] = df["hg_pval"].apply(lambda a: a if type(a) != str else np.max(np.array(a[1:-1].split(", "), dtype=np.float32)))

    filtered_go_ids=open(filtered_go_ids_file,'r').read().split() + [constants.ROOT_GO_ID]

    n_genes_pvals = [np.power([10 for a in range(x.count(",") + 1)], -np.array(x[1:-1].split(", ")).astype(np.float)) for x in df.reindex(filtered_go_ids).loc[:, "hg_pval"].dropna()]
    max_genes_pvals = reduce(lambda a, x: np.append(a, np.min(x)), n_genes_pvals, np.array([]))
    print("total n_genes with pval:{}/{}".format(max_genes_pvals.shape[0], len(filtered_go_ids)))
    max_genes_pvals = np.append(max_genes_pvals, np.ones(len(filtered_go_ids) - np.size(max_genes_pvals)))
    fdr_results = fdrcorrection0(max_genes_pvals, alpha=hg_th, method='indep', is_sorted=False)
    max_true_counter = np.sum(fdr_results[0])
    HG_CUTOFF = -np.log10(np.sort(max_genes_pvals))[max_true_counter - 1] if max_true_counter > 0 else -1
    print("HG cutoff: {}".format(HG_CUTOFF))

    included_terms=df.loc[filtered_go_ids, "hg_pval_max"].apply(lambda a: True)# apply(lambda row: np.any(np.array(row[1:-1].split(', ') if type(row) == str else [row],dtype=np.float) >= HG_CUTOFF))
    included_terms=df.index.isin(included_terms[included_terms].index)
    df_filtered_in = df.loc[included_terms, :]
    df_filtered_out = df.loc[~included_terms, :]

    if df_filtered_in.shape[0] != 0:
        df_filtered_in.loc[:,"hg_rank"] = df_filtered_in["hg_pval_max"].rank(ascending=0)
    else:
        df_filtered_in = pd.DataFrame(columns=['hg_rank'])

    df_all = pd.concat((df_filtered_in, df_filtered_out), axis=0)

    return df_all



def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--report_folder', dest='report_folder', default=constants.config_json["report_folder"])
    parser.add_argument('--filtered_go_ids_file', dest='filtered_go_ids_file', default=constants.config_json["filtered_go_ids_file"])
    parser.add_argument('--hg_th', dest='hg_th', default=constants.config_json["hg_th"])
    parser.add_argument('--n_permutations', dest='n_permutations', default=constants.config_json["n_permutations"])
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    report_folder=args.report_folder
    go_folder=args.go_folder
    filtered_go_ids_file=args.filtered_go_ids_file
    n_permutations=int(args.n_permutations)
    hg_th=float(args.hg_th)
    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    init_state(go_folder)

    print("{}_{}".format(dataset_name, algo))

    csv_file_name=os.path.join(report_folder, "emp_diff_modules_{dataset}_{algo}.tsv")

    df_all=add_md_to_terms(dataset_name, algo, n_permutations, csv_file_name=csv_file_name, filtered_go_ids_file=filtered_go_ids_file, hg_th=hg_th)
    try:
        df_all=df_all.loc[:, ["GO name", "hg_pval", "hg_pval_max", "hg_rank", "n_genes"]]
    except Exception:
        print("agg file was missing. No md file was produced")
        return 

    df_all.to_csv(csv_file_name.format(dataset=dataset_name, algo=algo)[:-4] + "_md.tsv", sep='\t') # df_all["hg_pval_max"].values > 0 # "passed_fdr", "emp_pval", "emp_rank"


if __name__=="__main__":
    main()

