import sys
sys.path.insert(0, '../')

import os
import subprocess
import src.constants as constants

import pandas as pd
import argparse


current_path=os.path.dirname(os.path.realpath(__file__))


def execute_stage(py_script, params):
    try:
        df_status_report = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t", index_col=0)
    except:
        print("no status report found. creating new one...")
        df_status_report=pd.DataFrame()

    df_status_report.loc[algo, dataset] = py_script
    df_status_report.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t")

    params = " ".join(params)
    print("about to start script {} with params:\n{}".format(py_script, params))
    out=subprocess.Popen("{}/../python27/bin/python {} {}".format(constants.REPO_DIR, py_script, params), shell=True,
                           stdout=subprocess.PIPE, cwd=current_path)
    print(out.stdout.read())
    out.wait()
    if out.returncode!=0:
        raise (Exception, "Error in {}: expected return code 0 but got {}".format(py_script, out.returncode))

    return out.returncode



def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset', dest='dataset', default="TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM")
    parser.add_argument('--omic_type', dest='omic_type', default="GE")
    parser.add_argument('--algo', dest='algo', default="my_netbox_td")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=5100)
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=20)
    parser.add_argument('--calc_true_scores', dest='calc_true_scores', default="false")
    parser.add_argument('--n_iteration', dest='n_iteration', default=1)
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=5000)
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=5000)
    parser.add_argument('--override_permutations', help="override_permutations", dest='override_permutations', default='false')
    parser.add_argument('--processes', help="processes", dest='aggregate_bg_distribution,add_go_metadata,calculate_significance', default=None)
    args = parser.parse_args()

    dataset=args.dataset.split(",")
    algo=args.algo.split(",")
    omic_type = args.omic_type
    network_file_name = args.network
    n_start=int(args.n_start)
    n_end=int(args.n_end)
    calc_true_scores=args.calc_true_scores.lower()
    override_permutations=args.override_permutations
    n_iteration = int(args.n_iteration)
    n_total_samples = int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    n_permutations = int(args.n_total_samples)
    pf=args.pf
    processes = args.processes

    omic_type_param = "--omic_type {}".format(omic_type)
    n_start_param = "--n_start {}".format(n_start)
    n_end_param = "--n_end {}".format(n_end)
    pf_param = "--pf {}".format(pf)
    calc_true_scores_param = "--calc_true_scores {}".format(calc_true_scores)
    n_iteration_param = "--n_iteration {}".format(n_iteration)
    n_total_samples_param = "--n_total_samples {}".format(n_total_samples)
    n_dist_samples_param = "--n_dist_samples {}".format(n_dist_samples)
    override_permutations_param = "--override_permutations {}".format(override_permutations)
    n_permutations_param = "--n_permutations {}".format(n_permutations)

    dataset_param="--dataset {}".format(dataset)
    algo_param="--algo {}".format(algo)

    params_by_processes={"generate_permuted_solutions": [dataset_param, algo_param, omic_type_param, n_start_param, n_end_param, pf_param, override_permutations_param],
                      "aggregate_bg_distribution": [dataset_param, algo_param, omic_type_param, n_start_param, n_end_param, pf_param, calc_true_scores_param],
                      "add_go_metadata": [dataset_param, algo_param, omic_type_param, n_permutations_param],
                      "calculate_significance": [dataset_param, algo_param, omic_type_param, pf_param, n_iteration_param, n_total_samples_param, n_dist_samples_param, n_iteration_param]}

    for cur_process in processes:
        if cur_process not in params_by_processes:
            print("unknown process detected: {}. abort...".format(cur_process))
            raise Exception

    try:
       for cur_process in processes:
           execute_stage(cur_process+".py", params_by_processes[cur_process])

    except Exception as e:
       print("error in {}, {}: {}".format(dataset, algo, e))
