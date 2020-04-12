import sys
sys.path.insert(0, '../..')

import src.constants as constants

import os
import numpy as np
import argparse

from src.utils.randomize_data_robustness import create_random_ds
from src.utils.randomize_data_robustness import permutation_dataset_exists
from src.utils.daemon_multiprocessing import MyPool, func_star


def empirical_dist_iteration(dataset_file, rand_idx, algo, network_file_name, ss_ratio, output_folder):

    # print("starting generate permuted datasets: {}, {}@".format(dataset_file, algo, rand_idx))
    create_random_ds(output_folder, dataset_file, network_file_name, ss_ratio, rand_idx)
    print("done generating permuted dataset: {}, {}, {}, {}@".format(dataset_file, algo, rand_idx, ss_ratio))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', default=constants.config_json["network_file"])
    parser.add_argument('--ss_ratios', dest='ss_ratios', default=constants.config_json["ss_ratios"])
    parser.add_argument('--robustness_datasets_folder', dest='robustness_datasets_folder', default=constants.config_json["robustness_datasets_folder"])
    parser.add_argument('--n_start_r', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start_r', default=constants.config_json["n_start_r"])
    parser.add_argument('--n_end_r', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end_r', default=constants.config_json["n_end_r"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--override_permutations', help="override_permutations", dest='override_permutations', default=constants.config_json["override_permutations"])

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    permuted_datasets_folder = args.robustness_datasets_folder

    ss_ratios = [float(a) for a in args.ss_ratios.split(",")]
    network_file = args.network_file

    parallelization_factor = int(args.pf)
    n_start=args.n_start_r
    n_end=args.n_end_r
    override_permutations=args.override_permutations.lower()=="true"

    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    break_loop=False
    while not break_loop:
        try:
            # [ empirical_dist_iteration(dataset_file, x, algo, network_file, ss_ratio, permuted_datasets_folder) for ss_ratio in ss_ratios for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_dataset_exists(dataset_name, x, ss_ratio, permuted_datasets_folder)]
            # empirical_dist_iteration(dataset_file, x, algo, output_folder)
            p = MyPool(parallelization_factor)
            params=[ [empirical_dist_iteration, [dataset_file, x, algo, network_file, ss_ratio, permuted_datasets_folder]] for ss_ratio in ss_ratios for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_dataset_exists(dataset_name, x, ss_ratio, permuted_datasets_folder)]
            p.map(func_star, params)
            break_loop=True

        except (MemoryError, OSError) as e:
            # p.close()
            print(e)
            raise
            # pass


if __name__ == "__main__":
    main()
           


