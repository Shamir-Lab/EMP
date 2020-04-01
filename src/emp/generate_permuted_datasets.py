import sys
sys.path.insert(0, '../')

import os
import numpy as np
import argparse

from src.utils.randomize_data import create_random_ds
from src.utils.randomize_data import permutation_dataset_exists
from src.utils.daemon_multiprocessing import MyPool, func_star


def empirical_dist_iteration(dataset_file, rand_idx, algo, output_folder):

    print("starting generate permuted datasets: {}, {}".format(dataset_file, algo, rand_idx))
    create_random_ds(output_folder, dataset_file, rand_idx)
    print("done generating permuted dataset: {}, {}, {}".format(dataset_file, algo, rand_idx))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data1/emp_test/datasets/brca.tsv")
    parser.add_argument('--algo', dest='algo', default="DOMINO")
    parser.add_argument('--permuted_datasets_folder', dest='permuted_datasets_folder', default="/media/hag007/Data/emp_test/output")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=5)
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=3)
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations', default="false")

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    permuted_datasets_folder = args.permuted_datasets_folder

    parallelization_factor = int(args.pf)
    n_start=args.n_start
    n_end=args.n_end
    override_permutations=args.override_permutations.lower()=="true"

    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    break_loop=False
    while not break_loop:
        try:
            # [ empirical_dist_iteration(dataset_file, x, algo, permuted_datasets_folder) for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_output_exists(dataset_file, algo, x, permuted_datasets_folder)]
            # empirical_dist_iteration(dataset_file, x, algo, output_folder)
            p = MyPool(parallelization_factor)
            params=[ [empirical_dist_iteration, [dataset_file, x, algo, permuted_datasets_folder]] for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_dataset_exists(dataset_name, x, permuted_datasets_folder)]
            p.map(func_star, params)
            break_loop=True

        except (MemoryError, OSError) as e:
            # p.close()
            print(e)
            raise
            # pass


if __name__ == "__main__":
    main()
           


