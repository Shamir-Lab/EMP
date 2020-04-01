import sys
sys.path.insert(0, '../')

import os
import numpy as np
import argparse
import json

from src.utils.randomize_data import get_permuted_folder_name
from src.utils.randomize_data import permutation_solution_exists
from src.utils.daemon_multiprocessing import MyPool, func_star
from src.runners.run_algo import run_algo
import src.constants as constants

def empirical_dist_iteration(dataset_file, rand_idx, algo, network_file, go_folder, permuted_datasets_folder, permuted_solutions_folder,  additional_args):

    print("starting generate permuted solution: {}, {}, {}".format(dataset_file, algo, rand_idx))
    dataset=os.path.splitext(os.path.split(dataset_file)[1])[0]
    permuted_folder=get_permuted_folder_name(dataset, index=rand_idx)
    output_folder=os.path.join(permuted_solutions_folder, "sol_{}_{}".format(algo,permuted_folder))
    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    run_algo(os.path.join(os.path.join(permuted_datasets_folder, permuted_folder), "data", constants.SCORES_FILE_NANE), algo=algo,
             network_file_name=network_file, go_folder=go_folder, output_folder=output_folder, **additional_args)
    print("done generating permuted solution: {}, {}, {}".format(dataset_file, algo, rand_idx))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data/emp_test/datasets/brca.tsv")
    parser.add_argument('--algo', dest='algo', default="DOMINO")
    parser.add_argument('--network_file', dest='network_file', help='/path/to/network_file', default="/media/hag007/Data/emp_test/networks/dip.sif")
    parser.add_argument('--go_folder', dest='go_folder', default="/media/hag007/Data/emp_test/go")
    parser.add_argument('--permuted_datasets_folder', dest='permuted_datasets_folder', default="/media/hag007/Data/emp_test/permuted_datasets")
    parser.add_argument('--permuted_solutions_folder', dest='permuted_solutions_folder', default="/media/hag007/Data/emp_test/permuted_solutions")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=5)
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=3)
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations', default="false")
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default='{"slices_file": "/media/hag007/Data1/emp_test/networks/dip_ng_modularity_components.txt", "modules_threshold": 0.05, "slices_threshold" : 0.3}')
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    network_file = args.network_file
    permuted_datasets_folder = args.permuted_datasets_folder
    permuted_solutions_folder = args.permuted_solutions_folder
    go_folder = args.go_folder
    additional_args=json.loads(args.additional_args)

    parallelization_factor =  int(args.pf)
    n_start=args.n_start
    n_end=args.n_end
    override_permutations=args.override_permutations.lower()=="true"

    constants.GO_DIR=go_folder

    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    break_loop=False
    while not break_loop:
        try:

            # [empirical_dist_iteration(dataset_file, x, algo, network_file, go_folder, permuted_datasets_folder, permuted_solutions_folder, additional_args) for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_solution_exists(dataset_name, algo, x, permuted_solutions_folder)]

            p = MyPool(parallelization_factor)
            params=[ [empirical_dist_iteration, [dataset_file, x, algo, network_file, go_folder, permuted_datasets_folder, permuted_solutions_folder, additional_args]] for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_solution_exists(dataset_name, algo, x, permuted_solutions_folder)]
            p.map(func_star, params)
            break_loop=True

        except (MemoryError, OSError) as e:
            # p.close()
            print(e)
            raise
            # break_loop=True
            # pass



if __name__ == "__main__":
    main()
           


