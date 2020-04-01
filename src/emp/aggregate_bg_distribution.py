import sys
sys.path.insert(0, '../')

import argparse

import pandas as pd
import numpy as np
import multiprocessing
import os

from pandas.errors import EmptyDataError

from src.utils import goids2gonames
from src.utils.daemon_multiprocessing import func_star
from src.utils.randomize_data import get_permuted_folder_name


def get_best_module_sig_score(report_folder, shared_list):

    try:

        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"

        print("fetch permutation {}".format(report_folder))
        n_modules=len(pd.read_csv(os.path.join(report_folder, "modules_summary.tsv"), sep='\t').index)
        go_results = [os.path.join(report_folder, cur_module) for cur_module in os.listdir(report_folder) if
                      "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]

        for cur in go_results:
            try:
                df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)
                df_go_pvals[df_go_pvals.isna()] = 1
                df_go_pvals=df_go_pvals.min(axis=1).to_frame()

            except EmptyDataError as e:
                print(e)
                pass

        if len(df_go_pvals.index)==0:
            df_go_pvals=pd.DataFrame(data=np.array([[1]]),index=["GO:0008150"])


        if shared_list is not None:
            shared_list.append(df_go_pvals)
            print("done aggregate {} permutations".format(len(shared_list)))

        return df_go_pvals

    except Exception as e:
           print(Exception, e)
           pass


def get_all_modules_sig_scores(report_folder, shared_list=None):

    try:

        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"

        print("fetch permutation {}".format(report_folder))
        n_modules=len(pd.read_csv(os.path.join(report_folder, "modules_summary.tsv"), sep='\t').index)
        go_results = [os.path.join(report_folder, cur_module) for cur_module in os.listdir(report_folder) if
                      "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]

        for cur in go_results:
            try:
                df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)

            except EmptyDataError as e:
                print(e)
                pass

        if len(df_go_pvals.index)==0:
            df_go_pvals=pd.DataFrame(data=np.array([[1]]),index=["GO:0008150"])

        df_go_pvals[df_go_pvals.isna()] = 1

        if shared_list is not None:
            shared_list.append(df_go_pvals)
            print("done aggregate {} permutations".format(len(shared_list)))

        return df_go_pvals

    except Exception as e:
           print(Exception, e)
           pass


def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data/emp_test/datasets/brca.tsv")
    parser.add_argument('--algo', dest='algo', default="DOMINO")
    parser.add_argument('--permuted_solutions_folder', dest='permuted_solutions_folder', default="/media/hag007/Data/emp_test/output")
    parser.add_argument('--report_folder', dest='report_folder', default="/media/hag007/Data/emp_test/output")
    parser.add_argument('--go_folder', dest='go_folder', default="/media/hag007/Data1/emp_test/go")
    parser.add_argument('--true_solution_folder', dest='true_solution_folder', default="/media/hag007/Data/emp_test/true_solutions")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=5)
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=5)

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    permuted_solutions_folder=args.permuted_solutions_folder
    true_solution_folder = args.true_solution_folder
    report_folder=args.report_folder
    go_folder = args.go_folder
    n_start=int(args.n_start)
    n_end= int(args.n_end)
    pf=args.pf

    manager = multiprocessing.Manager()
    l_permutations_top_pvals = manager.list()
    p = multiprocessing.Pool(int(pf))

    params=[]
    for cur_idx in range(n_start, n_end):
        permuted_folder=get_permuted_folder_name(os.path.splitext(os.path.split(dataset_file)[1])[0], cur_idx)
        sol_output_folder = os.path.join(permuted_solutions_folder, "sol_{}_{}".format(algo,permuted_folder), "report")
        print (permuted_solutions_folder)
        params.append([get_best_module_sig_score, [sol_output_folder, l_permutations_top_pvals]])
        print(n_end)
    p.map(func_star, params)

    df_all_terms = pd.concat(list(l_permutations_top_pvals), axis=1)
    df_all_terms=df_all_terms.fillna(1)
    print("total # permutations: {}/{}".format(len(l_permutations_top_pvals), n_end-n_start))

    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
    df_real_pvals = get_all_modules_sig_scores(os.path.join(true_solution_folder, "{}_{}".format(dataset_name,algo), "report"))
    df_real_pvals=df_real_pvals.apply(lambda x: -np.log10(x))

    df_real_pvals_as_list = df_real_pvals.apply(lambda x: str(list(x)), axis=1).to_frame()
    df_real_pvals_as_list.index = df_real_pvals.index

    print("total # real terms: {}".format(len(df_real_pvals.index)))
    df_results=df_all_terms.apply(lambda row : str(list(-np.log10(row.values.astype(np.float)))), axis=1).to_frame()
    df_results.columns = ['dist_n_samples']
    missing_indices=set(df_real_pvals.index).difference(df_results.index)
    print("missing_indices: {}".format(len(missing_indices)))
    print("existing_indices: {}".format(len(df_results.index)))

    for idx in missing_indices:
        df_results.loc[idx, "dist_n_samples"]=str([0 for a in np.arange(n_start, n_end)])

    df_results['hg_pval']= df_real_pvals_as_list.iloc[:, 0]
    df_results["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_real_pvals.index), go_folder),
                                   index=df_real_pvals.index)
    df_results.to_csv(os.path.join(report_folder,"emp_diff_modules_{}_{}.tsv".format(dataset_name, algo)),  sep='\t', index_label="GO id")

    print("permutation shape: {}".format(df_all_terms))

if __name__ == "__main__":
    main()