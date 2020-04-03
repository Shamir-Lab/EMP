import sys
sys.path.insert(0, '../..')
import json
import argparse
import shutil
import os
from src.runners.run_algo import run_algo




def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default="/media/hag007/Data1/emp_test/original_datasets/brca.tsv")
    parser.add_argument('--algo', dest='algo', default="netbox")
    parser.add_argument('--network_file', dest='network_file', help='/path/to/network_file', default="/media/hag007/Data1/emp_test/networks/dip.sif")
    parser.add_argument('--go_folder', dest='go_folder', default="/media/hag007/Data1/emp_test/go")
    parser.add_argument('--true_solution_folder', dest='true_solution_folder', default="/media/hag007/Data1/emp_test/true_solutions/")
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default='{"slices_file": "/media/hag007/Data1/emp_test/networks/dip_ng_modularity_components.txt", "modules_threshold": 0.05, "slices_threshold" : 0.3}')

    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    network_file = args.network_file
    go_folder = args.go_folder
    true_solution_folder = args.true_solution_folder
    additional_args = args.additional_args

    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    output_folder=os.path.join(true_solution_folder, "{}_{}".format(dataset_name,algo))
    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    run_algo(dataset_file, algo, network_file, go_folder, os.path.join(output_folder) , **json.loads(additional_args))


if __name__ == "__main__":
    main()