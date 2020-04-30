import sys
sys.path.insert(0, '../../')
# sys.path.insert(0, '../')
# sys.path.insert(0, '../../../')

import os
import subprocess
import src.constants as constants

import pandas as pd
import argparse
import json

current_path=os.path.dirname(os.path.realpath(__file__))


def execute_stage(py_script, params):
    try:
        df_status_report = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t", index_col=0)
    except:
        pass
        # print("no status report found. creating new one...")
        df_status_report=pd.DataFrame()

    # df_status_report.loc[algo, dataset] = py_script
    # df_status_report.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t")

    params = " ".join(params)
    print("about to start script {} with params:\n{}".format(py_script, params))
    prc=subprocess.Popen("{}python/bin/python {} {}".format(constants.dir_path, py_script, params), shell=True,
                           stdout=subprocess.PIPE, cwd=current_path)
    while True:
        output = prc.stdout.readline()
        if output == b'' and prc.poll() is not None:
            break
        if output:
            out_str=output.decode("utf-8")
            if out_str.endswith("@\n"):
                print('\r'+out_str[:-2], end ='')
            else:
                print(out_str, end ='')

    rc = prc.poll()

    if rc!=0:
        raise Exception("Error in {}: expected return code 0 but got {}".format(py_script, rc))

    # os.path.join(constants.dir_path, "src", "emp")
    return rc



def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--permuted_datasets_folder', dest='permuted_datasets_folder', default=constants.config_json["permuted_datasets_folder"])
    parser.add_argument('--permuted_solutions_folder', dest='permuted_solutions_folder', default=constants.config_json["permuted_solutions_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--report_folder', dest='report_folder', default=constants.config_json["report_folder"])
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=constants.config_json["n_start"])
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=constants.config_json["n_end"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--n_permutations', dest='n_permutations', default=constants.config_json["n_permutations"])
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=constants.config_json["n_total_samples"])
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=constants.config_json["n_dist_samples"])
    parser.add_argument('--override_permutations', help="override_permutations", dest='override_permutations', default=constants.config_json["override_permutations"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default=constants.config_json["additional_args"])
    parser.add_argument('--processes', dest="processes", help='regenerate_bg_distribution,add_go_metadata,calculate_significance', default=constants.config_json["processes"])
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    network_file=args.network_file
    go_folder=args.go_folder
    permuted_datasets_folder=args.permuted_datasets_folder
    permuted_solutions_folder=args.permuted_solutions_folder
    true_solutions_folder=args.true_solutions_folder
    report_folder=args.report_folder
    n_start=int(args.n_start)
    n_end=int(args.n_end)
    pf=args.pf
    override_permutations=args.override_permutations
    n_permutations = int(args.n_total_samples)
    n_total_samples = int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    additional_args = args.additional_args
    processes = args.processes.split(",")

    dataset_file_params = "--dataset_file {}".format(dataset_file)
    algo_param="--algo {}".format(algo)
    permuted_datasets_folder_param = "--permuted_datasets_folder {}".format(permuted_datasets_folder)
    permuted_solutions_folder_param = "--permuted_solutions_folder {}".format(permuted_solutions_folder)
    true_solutions_folder_param = "--true_solutions_folder {}".format(true_solutions_folder)
    report_folder_param = "--report_folder {}".format(report_folder)
    network_file_param = "--network_file {}".format(network_file)
    go_folder_param = "--go_folder {}".format(go_folder)
    n_start_param = "--n_start {}".format(n_start)
    n_end_param = "--n_end {}".format(n_end)
    pf_param = "--pf {}".format(pf)
    override_permutations_param = "--override_permutations {}".format(override_permutations)
    n_permutations_param = "--n_permutations {}".format(n_permutations)
    n_total_samples_param = "--n_total_samples {}".format(n_total_samples)
    n_dist_samples_param = "--n_dist_samples {}".format(n_dist_samples)
    additional_args_param = "--additional_args {}".format(json.dumps(str(additional_args)))

    params_by_processes={
        "generate_permuted_datasets": [dataset_file_params, algo_param, permuted_datasets_folder_param, n_start_param, n_end_param, pf_param, override_permutations_param],
        "generate_permuted_solutions": [dataset_file_params, algo_param, network_file_param, go_folder_param, permuted_datasets_folder_param, permuted_solutions_folder_param, n_start_param, n_end_param, pf_param, override_permutations_param, additional_args_param],
        "generate_true_solution": [dataset_file_params, algo_param, network_file_param, go_folder_param, true_solutions_folder_param, additional_args_param],
        "aggregate_bg_distribution": [dataset_file_params, algo_param, go_folder_param, permuted_solutions_folder_param,  true_solutions_folder_param, report_folder_param, n_start_param, n_end_param, pf_param],
        "regenerate_bg_distribution": [dataset_file_params, algo_param, go_folder_param, permuted_solutions_folder_param,  true_solutions_folder_param, report_folder_param, n_start_param, n_end_param, pf_param],
        "add_go_metadata": [dataset_file_params, algo_param, go_folder_param, report_folder_param, n_permutations_param],
        "calculate_significance": [dataset_file_params, algo_param, report_folder_param, n_total_samples_param, n_dist_samples_param]}

    for cur_process in processes:
        if cur_process not in params_by_processes:
            print("unknown process detected: {}. abort...".format(cur_process))
            raise Exception

    try:
       for cur_process in processes:
           execute_stage(cur_process+".py", params_by_processes[cur_process])

    except Exception as e:
       print("error in {}: {}".format(cur_process, e))
       raise


if __name__=="__main__":
    main()

