import sys
sys.path.insert(0, '../')

import os

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from functools import reduce

from src import constants
from src.utils import go_hierarcies
from src.utils.ensembl2entrez import entrez2ensembl_convertor

# constants.GO_DIR="/media/hag007/Data/bnet/GO"
dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
    roots=['GO:0008150'])
vertices = list(dict_result.values())[0]['vertices']

terms_to_genes = {}


def mean_difference(row, dataset_data, classes_data):
    try:
        return dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                                True)), classes_data == 2].dropna().values.mean() - \
               dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                                True)), classes_data == 1].dropna().values.mean()

    except:
        print("no gene were found for {}, {} (pval={})".format(row["index"], row["GO name"],
                                                               row["hg_pval"]))


def calc_empirical_pval(row, n_permutation):
    emp_pvals = []
    if type(row["hg_pval"]) == str:
        hg_pvals = np.array(row["hg_pval"][1:-1].split(", "), dtype=np.float)
    else:
        hg_pvals = np.array([row['hg_pval']], dtype=np.float)

    for hg_pval in hg_pvals:

        pval = np.array([float(x) for x in row["dist_n_samples"][1:-1].split(", ")])

        if len(pval) < n_permutation:
            raise ValueError("too few samples: {} (expected at least {})".format(len(pval), n_permutation))

        else:
            hg_pval = np.array(hg_pval, dtype=np.float32)
            pval = np.array(pval[:n_permutation], dtype=np.float32)
            pos = np.size(pval) - np.searchsorted(np.sort(pval), hg_pval, side='left')
            if row['GO name'] == "lymphocyte migration":
                print("pos: {}, hg_pval: {}, max dist: {}".format(pos, hg_pval, np.max(pval)))
            emp_pval = pos / float(np.size(pval))
            emp_pvals.append(emp_pval)

    return str(emp_pvals)


def get_all_genes_for_term(vertices, cur_root, term, in_subtree):
    if term in terms_to_genes:
        return terms_to_genes[cur_root]

    all_genes = set()
    if in_subtree:
        try:
            all_genes.update(go2geneids[cur_root])
        except Exception:
            # print(e)
            pass
            # print("cur_root : {} was not found".format(cur_root))

    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))

    terms_to_genes[cur_root] = all_genes
    return all_genes


def add_md_to_terms(dataset="SOC", algo="jactivemodules_sa", n_permutations=300, csv_file_name=None):

    csv_file_name = csv_file_name.format(dataset=dataset, algo=algo)
    try:
        df = pd.read_csv(csv_file_name, sep='\t', index_col=0)
    except:
        print("no bg file {} for {}, {}".format(csv_file_name, dataset,algo))
        return None
    df = df.dropna()

    n_genes = [len(get_all_genes_for_term(vertices, cur_go_id, cur_go_id, cur_go_id == cur_go_id)) for i, cur_go_id in
               enumerate(df.index.values)]
    depth = [list(dict_result.values())[0]['vertices'][cur_go_id]['D'] for i, cur_go_id in enumerate(df.index.values)]
    df["n_genes"] = pd.Series(n_genes, index=df.index)
    df["depth"] = pd.Series(depth, index=df.index)
    df = df.rename(columns={"filtered_pval": "hg_pval"})
    df["hg_pval_max"] = df["hg_pval"].apply(
        lambda a: a if type(a) != str else np.max(np.array(a[1:-1].split(", "), dtype=np.float32)))

    n_genes_pvals = df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500]), "hg_pval"]
    n_modules = n_genes_pvals.shape[0]
    if n_modules > 0:
        n_modules = str(n_genes_pvals.iloc[0]).count(",") + 1
    n_genes_pvals = n_genes_pvals.values
    print("start reduction..")
    n_genes_pvals = [
        np.power([10 for a in range(x.count(",") + 1)], -np.array(x[1:-1].split(", ")).astype(np.float)) if type(
            x) == str else [10 ** (-x)] for x in n_genes_pvals]
    max_genes_pvals = reduce(lambda a, x: np.append(a, np.min(x)), n_genes_pvals, np.array([]))

    print("total n_genes with pval:{}/{}".format(max_genes_pvals.shape[0], constants.N_GO_TERMS))
    max_genes_pvals = np.append(max_genes_pvals, np.ones(constants.N_GO_TERMS - np.size(max_genes_pvals)))
    fdr_results = fdrcorrection0(max_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    max_true_counter = np.sum(fdr_results[0])
    HG_CUTOFF = -np.log10(np.sort(max_genes_pvals))[max_true_counter - 1] if max_true_counter > 0 else 0
    print("HG cutoff: {}".format(HG_CUTOFF))

    df_filtered_in = df.loc[np.logical_and.reduce([df.loc[:, "n_genes"].values > 5, df.loc[:, "n_genes"].values < 500,
                                                   df.loc[:, "hg_pval"].apply(lambda row: np.any(
                                                       np.array(row[1:-1].split(', ') if type(row) == str else [row],
                                                                dtype=np.float) >= HG_CUTOFF))]), :]
    df_filtered_out = df.loc[~np.logical_and.reduce([df.loc[:, "n_genes"].values > 5, df.loc[:, "n_genes"].values < 500,
                                                     df.loc[:, "hg_pval"].apply(lambda row: np.any(
                                                         np.array(row[1:-1].split(', ') if type(row) == str else [row],
                                                                  dtype=np.float) >= HG_CUTOFF))]), :]

    # if df_filtered_in.shape[0] != 0:
    #
    #     df_filtered_in.loc[:,"index"] = df_filtered_in.index.values
    #     df_filtered_in.loc[:,"emp_pval"] = df_filtered_in.apply(lambda row: calc_empirical_pval(row, n_permutations), axis=1)
    #     df_filtered_in.loc[:,"emp_pval_max"] = df_filtered_in.loc[:,"emp_pval"].apply(
    #         lambda a: a if type(a) != str else np.min(np.array(a[1:-1].split(", "), dtype=np.float32)))
    #
    #     pvals_corrected = [[x] if type(x) != str else np.array(x[1:-1].split(", ")).astype(np.float) for x in
    #                        df_filtered_in.loc[:, "emp_pval"]]
    #
    #     max_pvals_corrected = reduce(lambda a, x: np.append(a, np.min(x)), pvals_corrected, np.array([]))
    #
    #     print("# of corrected pvals: {}".format(max_pvals_corrected.shape[0]))
    #
    #     pvals_corrected = reduce(lambda a, x: np.append(a, x), pvals_corrected, np.array([]))
    #     fdr_results = fdrcorrection0(max_pvals_corrected, alpha=0.05, method='indep', is_sorted=False)
    #     max_true_counter = len([cur for cur in fdr_results[0] if cur == True])
    #     emp_cutoff = np.sort(np.sort(max_pvals_corrected))[max_true_counter - 1] if max_true_counter > 0 else 0
    #     print("number of true hypothesis: {} (emp cutoff: {}, n={})".format(max_true_counter, emp_cutoff, len(fdr_results[0])))
    #     if n_modules > 1:
    #         df_filtered_in.loc[:,"passed_fdr"] = [str(a) for a in (pvals_corrected <= emp_cutoff).reshape(-1, n_modules)]
    #
    #     else:
    #         df_filtered_in.loc[:,"passed_fdr"] = pvals_corrected <= emp_cutoff
    #         df_filtered_in.loc[:,"passed_fdr"] = df_filtered_in["passed_fdr"].apply(lambda a: str([a]))
    #
    #     df_filtered_in.loc[:,"emp_rank"] = df_filtered_in["emp_pval_max"].rank(ascending=1)
    #     df_filtered_in.loc[:,"hg_rank"] = df_filtered_in["hg_pval_max"].rank(ascending=0)
    #
    #     df_filtered_in = df_filtered_in.sort_values(by=["emp_rank", "hg_rank"])
    #
    # else:
    #     df_filtered_in = pd.DataFrame(columns=['emp_pval', 'hg_rank', 'emp_rank', 'passed_fdr'])

    if df_filtered_in.shape[0] != 0:
        df_filtered_in.loc[:,"hg_rank"] = df_filtered_in["hg_pval_max"].rank(ascending=0)

    else:
        df_filtered_in = pd.DataFrame(columns=['hg_rank']) # 'emp_pval', 'hg_rank', 'emp_rank', 'passed_fdr'

    df_all = pd.concat((df_filtered_in, df_filtered_out), axis=0)

    return df_all

# if __name__=='__main__':
#     get_all_genes_for_term(vertices, "GO:0008150", "GO:0008150", True)
#     x=1
