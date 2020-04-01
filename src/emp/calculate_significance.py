import sys
sys.path.insert(0, '../')

from functools import reduce
import argparse
from src import constants
import os
import pandas as pd
import numpy as np
import random
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import multiprocessing


def calc_emp_pval(cur_rv, cur_dist):
    cur_dist=np.array(cur_dist, np.float32)
    emp_pvals=[]
    if type(cur_rv)==str:
        hg_pvals=np.array(cur_rv[1:-1].split(", "), dtype=np.float32)
    else:
        hg_pvals=np.array([cur_rv], dtype=np.float32)

    for hg_pval in hg_pvals:
        
        pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), hg_pval, side='left')
        emp_pval = pos / float(np.size(cur_dist))
        emp_pvals.append(emp_pval)

    return str(emp_pvals)


def calculate_sig(algo_sample = None, dataset_sample = None, n_dist_samples = 300, n_total_samples = None, n_start_i=None, limit = 10000, md_path=None, dist_path=None):


    output_md = pd.read_csv(md_path.format(dataset_sample, algo_sample), sep='\t', index_col=0).dropna()
    output_md = output_md.rename(columns={"filtered_pval": "hg_pval"})
    n_genes_pvals=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), "hg_pval"]
    

    if n_genes_pvals.shape[0]==0:
       n_modules=0
       n_terms=0
       n_genes_pvals=np.array([])
    else:
       n_modules=str(n_genes_pvals.iloc[0]).count(",")+1
       n_terms=np.size(n_genes_pvals)
       n_genes_pvals=n_genes_pvals.values
    print("start reduction..")
    n_genes_pvals = [np.power([10 for a in range(x.count(",")+1)],-np.array(x[1:-1].split(", ")).astype(np.float)) if type(x)==str else [10**(-x)] for x in n_genes_pvals]
    max_genes_pvals = reduce(lambda a, x : np.append(a,np.min(x)), n_genes_pvals , np.array([]))


    print("total n_genes with pval:{}/{}".format(np.size(max_genes_pvals), 7435))
    max_genes_pvals=np.append(max_genes_pvals,np.ones(7435-np.size(max_genes_pvals)))
    fdr_results = fdrcorrection0(max_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=np.sort(max_genes_pvals)[n_hg_true-1]
    print("HG cutoff: {}, (n={})".format(HG_CUTOFF, n_hg_true))

    output_md = output_md.loc[np.logical_and.reduce(
        [output_md["n_genes"].values > 5, output_md["n_genes"].values < 500,
         output_md["hg_pval_max"].values > np.log10(HG_CUTOFF)]), :]

    output = pd.read_csv(dist_path.format(dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()
    output = output.rename(columns={"filtered_pval": "hg_pval"})
    output = output.loc[output_md.index.values, :]
    counter = 0
    emp_dists = []
    emp_pvals = []

    n_total_samples=n_total_samples if n_total_samples is not None else len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
    if n_start_i is None:
        np.random.seed(int(random.random()*1000))        
        i_choice=np.random.choice(n_total_samples, n_dist_samples, replace=False)
        i_dist=i_choice[:n_dist_samples]
    else:
        i_dist=np.arange(n_start_i, n_start_i+n_dist_samples)



    for index, cur in output.iterrows():
        if counter == limit: break
        pval = np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])[i_dist]
            
        emp_pvals.append(calc_emp_pval(cur["hg_pval"], pval))
        output_md.loc[index,'emp_pval']=emp_pvals[-1]
        emp_dists.append(pval)
        counter += 1

    mask_ids=output.index.values    
    emp_pvals = [np.array([x]) if type(x)!=str else np.array(x[1:-1].split(", ")).astype(np.float32)
                       for x in emp_pvals]
 
    if len(emp_pvals)==0:
        n_modules=0
        n_terms=0
    else:
        n_modules=emp_pvals[0].shape[0]
        n_terms=len(emp_pvals)

    max_emp_pvals = reduce(lambda a, x: np.append(a, np.min(x)), emp_pvals, np.array([]))
    emp_pvals = reduce(lambda a, x: np.append(a, x), emp_pvals, np.array([]))
    dist_name = "emp"
    df_dists = pd.DataFrame(index=output.index)
    df_dists["emp"] = pd.Series(emp_dists, index=output.index[:limit])
    max_emp_pvals = np.sort([x if x != 0 else 1.0 / n_dist_samples for x in max_emp_pvals])
    print("max emp pvals len: {}".format(len(max_emp_pvals)))
    print("min vals", 1.0 / n_dist_samples, np.min(list(max_emp_pvals) + [1]))
    print("max_genes_pvals: {}".format(max_genes_pvals.shape[0]))
    fdr_bh_results = fdrcorrection0(max_emp_pvals, alpha=0.05, method='indep', is_sorted=False)[0]
    n_emp_true=np.sum(fdr_bh_results)
    print("n_emp_true: {}".format(n_emp_true))

    if n_emp_true==0: 
       EMP_TH=0
    else: 
       EMP_TH=np.sort(max_emp_pvals)[n_emp_true-1]
     
    mask_terms=np.array([max(a, 1.0 / n_dist_samples)<=EMP_TH for a in emp_pvals])
    go_ids_result=np.array([])
    go_names_result=np.array([])
    n_emp_true=0
    if len(mask_terms) > 0:
       mask_terms=np.array(mask_terms).reshape(-1, n_modules)
       go_ids_result=output.index.values[mask_terms.any(axis=1)]
       go_names_result=output["GO name"].values[mask_terms.any(axis=1)]
       n_emp_true =np.sum(mask_terms)


    print("EMP cutoff: {} # true terms passed EMP cutoff: {}".format(EMP_TH, n_emp_true))

    return EMP_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms


def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data1/emp_test/datasets/scz.tsv")
    parser.add_argument('--algo', dest='algo', default="DOMINO")
    parser.add_argument('--report_folder', dest='report_folder', default="/media/hag007/Data1/emp_test/report")
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=10)
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=10)


    args = parser.parse_args()
    dataset_file=args.dataset_file
    algo=args.algo
    n_total_samples=int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    report_folder = args.report_folder

    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    print("about to calculate pval for {},_{}".format(dataset_name, algo))


    md_path=os.path.join(report_folder, "emp_diff_modules_{}_{}_md.tsv")
    dist_path=os.path.join(report_folder,"emp_diff_modules_{}_{}.tsv")

    BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result, mask_ids, mask_terms = \
        calculate_sig(algo, dataset_name, n_dist_samples, n_total_samples, md_path=md_path, dist_path=dist_path)


    output_md = pd.read_csv(
        os.path.join(report_folder, "emp_diff_modules_{}_{}_md.tsv".format(dataset_name, algo)),
        sep='\t', index_col=0)

    output_md.loc[mask_ids,"passed_oob_permutation_test"]="[False]"
    output_md.loc[mask_ids,"passed_oob_permutation_test"]=[str(list(a)) for a in mask_terms]
    print("emp_pval:\n", output_md['emp_pval'].iloc[0], type(output_md['emp_pval'].iloc[0]))

    if type(output_md['emp_pval'].iloc[0])!=str:
        output_md['emp_pval_max']=1.0
    else:
        print("output_md :{}".format(np.min(np.array(output_md['emp_pval'].iloc[0][1:-1].split(", "),dtype=np.float))))
        output_md['emp_pval_max']=output_md['emp_pval'].apply(lambda a: np.min(np.array(a[1:-1].split(", "),dtype=np.float)) if type(a)==str else "" )

    output_md.to_csv(
        os.path.join(report_folder, "emp_diff_modules_{}_{}_passed_oob.tsv".format(dataset_name, algo)),
        sep='\t')


if __name__ == "__main__":
    main()

