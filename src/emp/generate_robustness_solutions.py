import sys
sys.path.insert(0, '../..')

import os
import numpy as np
import argparse
import json

from src.utils.randomize_data_robustness import get_permuted_folder_name
from src.utils.randomize_data_robustness import permutation_solution_exists
from src.utils.daemon_multiprocessing import MyPool, func_star
from src.runners.run_algo import run_algo
import src.constants as constants
from src.utils.go import init_state

def empirical_dist_iteration(dataset_file, rand_idx, algo, ss_ratio, network_file, go_folder, permuted_datasets_folder, permuted_solutions_folder,  additional_args):

    print("starting generate permuted solution: {}, {}, {}, {}".format(dataset_file, algo, rand_idx, ss_ratio))
    dataset=os.path.splitext(os.path.split(dataset_file)[1])[0]
    permuted_folder=get_permuted_folder_name(dataset, index=rand_idx, ss_ratio=ss_ratio)
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
    parser.add_argument('--dataset_file', dest='dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--robustness_datasets_folder', dest='robustness_datasets_folder', default=constants.config_json["robustness_datasets_folder"])
    parser.add_argument('--robustness_solutions_folder', dest='robustness_solutions_folder', default=constants.config_json["robustness_solutions_folder"])
    parser.add_argument('--n_start_r', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start_r', default=constants.config_json["n_start_r"])
    parser.add_argument('--n_end_r', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end_r', default=constants.config_json["n_end_r"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--ss_ratios', help="ss_ratios", dest='ss_ratios', default=constants.config_json["ss_ratios"])
    parser.add_argument('--override_permutations', help="override_permutations", dest='override_permutations', default=constants.config_json["override_permutations"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default="{\"slices_file\": \"/media/hag007/Data/emp_test/networks/dip_ng_modularity_components.txt\", \"modules_threshold\": 0.05, \"slices_threshold\" : 0.3}")

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    ss_ratios = [float(a) for a in args.ss_ratios.split(",")]
    network_file = args.network_file
    robustness_datasets_folder = args.robustness_datasets_folder
    robustness_solutions_folder = args.robustness_solutions_folder
    go_folder = args.go_folder
    additional_args=json.loads(args.additional_args)
    parallelization_factor =  int(args.pf)
    n_start=args.n_start_r
    n_end=args.n_end_r
    override_permutations=args.override_permutations.lower()=="true"
    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]

    init_state(go_folder)

    break_loop=False
    while not break_loop:
        try:

            # [empirical_dist_iteration(dataset_file, x, algo, network_file, go_folder, permuted_datasets_folder, permuted_solutions_folder, additional_args) for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_solution_exists(dataset_name, algo, x, permuted_solutions_folder)]

            p = MyPool(parallelization_factor)
            params=[ [empirical_dist_iteration, [dataset_file, x, algo, ss_ratio, network_file, go_folder, robustness_datasets_folder, robustness_solutions_folder, additional_args]] for ss_ratio in ss_ratios for x in np.arange(int(n_start), int(n_end)) if override_permutations or not permutation_solution_exists(dataset_name, algo, x, ss_ratio, robustness_solutions_folder)]
            print("about to start generation of {} permuted solutions".format(len(params)))
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


