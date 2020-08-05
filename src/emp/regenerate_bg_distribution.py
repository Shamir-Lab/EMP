import sys
sys.path.insert(0, '../..')

import argparse

import pandas as pd
import numpy as np
import multiprocessing
import os
import src.constants as constants

from pandas.errors import EmptyDataError

from src.utils import go
from src.utils.daemon_multiprocessing import func_star
from src.utils.randomize_data import get_permuted_folder_name
from src.utils.go import init_state

def get_best_module_sig_score(report_folder, shared_list):

    try:

        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"

        # print("@fetch permutation {}".format(report_folder))
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
            df_go_pvals=pd.DataFrame(data=np.array([[1]]), index=[constants.ROOT_GO_ID])


        if shared_list is not None:
            shared_list.append(df_go_pvals)
            print("done aggregate {} permutations".format(len(shared_list)), end="@\n")

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
    parser.add_argument('--dataset_file', dest='dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--permuted_solutions_folder', dest='permuted_solutions_folder', default=constants.config_json["permuted_solutions_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--report_folder', dest='report_folder', default=constants.config_json["report_folder"])
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=constants.config_json["n_start"])
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=constants.config_json["n_end"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    permuted_solutions_folder=os.path.abspath(args.permuted_solutions_folder)
    true_solutions_folder = os.path.abspath(args.true_solutions_folder)
    report_folder=args.report_folder
    go_folder = args.go_folder
    n_start=int(args.n_start)
    n_end= int(args.n_end)
    pf=args.pf

    init_state(go_folder)


    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    df_results=pd.read_csv(os.path.join(report_folder,"emp_diff_modules_{}_{}_original.tsv".format(dataset_name, algo)),  sep='\t', index_col=0)
    df_results=df_results.loc[:,["dist_n_samples"]]


    df_real_pvals = get_all_modules_sig_scores(os.path.join(true_solutions_folder, "{}_{}".format(dataset_name,algo), "report"))
    df_real_pvals=df_real_pvals.apply(lambda x: -np.log10(x))

    df_real_pvals_as_list = df_real_pvals.apply(lambda x: str(list(x)), axis=1).to_frame()
    df_real_pvals_as_list.index = df_real_pvals.index


    df_results['hg_pval']= df_real_pvals_as_list.iloc[:, 0]
    df_results["GO name"] = pd.Series(go.get_go_names(list(df_real_pvals.index), go_folder),
                                   index=df_real_pvals.index)
    df_results.to_csv(os.path.join(report_folder,"emp_diff_modules_{}_{}.tsv".format(dataset_name, algo)),  sep='\t', index_label="GO id")


if __name__ == "__main__":
    main()
