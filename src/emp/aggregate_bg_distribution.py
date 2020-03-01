import sys
sys.path.insert(0, '../')

import argparse

import pandas as pd
import numpy as np
import shutil
import multiprocessing
import os

from pandas.errors import EmptyDataError

from src import constants
from src.utils import goids2gonames
from src.utils.daemon_multiprocessing import func_star
from src.utils.randomize_data import get_permutation_name
from src.runners.datasets_multithread_runner import run_dataset

def calc_dist(algos, datasets, shared_list=None, is_max=True):
    try:

        for cur_algo in algos:
            algos_filter = cur_algo

            df_go_pvals = pd.DataFrame()
            df_go_pvals.index.name="GO id"
            for cur_ds in datasets:
                print("fetch permutation {}".format(cur_ds))
                n_modules=len(pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "modules_summary.tsv"), sep='\t').index)
                go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                              os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                              if os.path.isdir(
                        os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo ==  algos_filter for
                              cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                              "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]

                for cur in go_results:
                    try:
                        df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)

                        if is_max:
                            df_go_pvals[df_go_pvals.isna()] = 1
                            df_go_pvals=df_go_pvals.min(axis=1).to_frame()

                    except EmptyDataError as e:
                        print(e)
                        pass
                if len(go_results)==0:
                    df_go_pvals=pd.DataFrame(data=np.array([[1]]),index=["GO:0008150"])

            if not is_max:
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
    parser.add_argument('--dataset', dest='dataset', default="SOC")
    parser.add_argument('--algo', dest='algo', default="jactivemodules_greedy")
    parser.add_argument('--omic_type', dest='omic_type', default="GE")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--output_folder', dest='output_folder', default=None)
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=100)
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=10)
    parser.add_argument('--calc_true_scores', dest='calc_true_scores', default="true")

    args = parser.parse_args()

    dataset=args.datasets
    algo=args.algos
    omic_type = args.omic_type
    network_file_name = args.network
    output_folder = args.output_folder
    if output_folder is None:
        output_folder=os.path.join(dataset,"output")
    n_start=int(args.n_start)
    n_end=int(args.n_end)
    calc_true_scores=args.calc_true_scores.lower()=="true"
    pf=args.pf


    for dataset in dataset:

        score_method = constants.PREDEFINED_SCORE

        for algo in algo:
            manager = multiprocessing.Manager()
            pvals = manager.list()
            p = multiprocessing.Pool(int(pf))
            params=[[calc_dist, [[algo], [get_permutation_name(omic_type, dataset, algo,  cur)], pvals]] for cur in range(n_start,n_end)]

            p.map(func_star, params)
            pvals=list(pvals)
            df_all_terms = pd.concat(pvals, axis=1)
            df_all_terms=df_all_terms.fillna(1)
            print("total # permutations: {}/{}".format(len(pvals), n_end-n_start))
            
            if calc_true_scores:
                if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,  "{}_{}".format(omic_type, dataset), algo)):
                    shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(omic_type, dataset), algo))
                run_dataset("{}_{}".format(omic_type, dataset), score_method=score_method,
                            algos=[algo], network_file_name=network_file_name)

            df_real_agg_pval=df_real_agg_pval.apply(lambda x: -np.log10(x))
            df_real_max_pval=df_real_agg_pval.max(axis=1).to_frame()

            df_real_agg_list_pval = df_real_agg_pval.apply(lambda x: str(list(-np.log10(x))), axis=1).to_frame()
            df_real_agg_list_pval.index = df_real_agg_pval.index

            print("total # real terms: {}".format(len(df_real_max_pval.index)))
            df_results=df_all_terms.apply(lambda row : str(list(-np.log10(row.values.astype(np.float)))), axis=1).to_frame()
            df_results.columns = ['dist_n_samples']
            missing_indices=set(df_real_max_pval.index).difference(df_results.index)
            df_real_agg_list_pval.loc[missing_indices, "dist_n_samples"]=str([0])
            df_results['hg_pval']= df_real_agg_list_pval.iloc[:, 0]
            df_results["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_real_max_pval.index)),
                                           index=df_real_max_pval.index)
            df_results.loc[df_results["hg_pval"].isna().values, "hg_pval"] = 0

            df_results.to_csv(os.path.join(output_folder,"emp_diff_modules_{}_{}.tsv".format(dataset, algo)),  sep='\t', index_label="GO id")

            print("permutation shape: {}".format(df_all_terms))
   

if __name__ == "__main__":
    main()