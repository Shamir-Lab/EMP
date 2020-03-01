import sys
sys.path.insert(0, '../')

import os
import numpy as np
import argparse
import pandas as pd
from pandas.errors import EmptyDataError

from src import constants
from src.utils.randomize_data import create_random_ds
from src.utils.randomize_data import permutation_output_exists
from src.utils.daemon_multiprocessing import MyPool, func_star
from src.runners.datasets_multithread_runner import run_dataset

def calc_dist(algos, datasets,empirical_th=None):
    for cur_algo in algos:
        algos_filter = cur_algo

        df_go = pd.DataFrame(columns=['qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                          os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for
                          cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                          "separated_modules" in cur_module]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                    df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)
                except EmptyDataError:
                    pass
        df_go_pvals[df_go_pvals.isna()]=1
        df_go = df_go[df_go['qval'] < 0.05]
        if empirical_th:
            df_go = df_go[df_go['pval'].apply(lambda x:-np.log10(x)) > empirical_th]

        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        return pval, df_go, df_go_pvals


def empirical_dist_iteration(prefix, dataset, cur, algo, score_method, output_folder, network_file_name="dip.sif"):

    print("starting iteration: {}, {}, {}".format(prefix, dataset, cur))
    random_ds = create_random_ds(prefix, "{}_{}".format(prefix, dataset), cur, algo)
    permuted_network_file_name = network_file_name
    run_dataset(random_ds, score_method=score_method,
                algos=[algo], network_file_name=permuted_network_file_name, output_folder=output_folder)
    cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
    print("done iteration: {}, {}, {}".format(prefix, dataset, cur))
    return cur_pval, df_pval_terms



def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset', dest='dataset', default="SOC")
    parser.add_argument('--omic_type', dest='omic_type', default="PASCAL_SUM")
    parser.add_argument('--algo', dest='algo', default="jactivemodules_greedy")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--output_folder', dest='output_folder', default=None)
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=100)
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=10)
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations', default="false")

    args = parser.parse_args()

    dataset=args.dataset
    algo=args.algo
    omic_type = args.omic_type
    network_file_name = args.network
    output_folder = args.output_folder
    if output_folder is None:
        output_folder=os.path.join(os.path.split(dataset)[0],"output")

    parallelization_factor = int(args.pf)
    n_start=args.n_start
    n_end=args.n_end
    override_permutations=args.override_permutations.lower()=="true"

    for dataset in dataset:

        score_method = constants.PREDEFINED_SCORE

        for algo in algo:
            break_loop=False
            while not break_loop:
                try:
                    p = MyPool(parallelization_factor)
                    params=[ [empirical_dist_iteration, [omic_type, dataset, x, algo, score_method, output_folder, network_file_name]] for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_output_exists(omic_type, dataset, algo, x)]
                    p.map(func_star, params)
                    break_loop=True

                except (MemoryError, OSError) as e:
                    p.close()
                    pass



if __name__ == "__main__":
    main()
           


