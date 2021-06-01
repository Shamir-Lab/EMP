# EMP

EMP: EMpirical Pipeline for correcting GO terms obtained by active module identification (AMI) algorithms.

AMI algorithms receive a gene network and nodes' activity scores as input and report sub-networks (modules) that are putatively biologically active. We observed that GO terms enriched in modules detected by these methods on the real data were often also enriched after randomly permuting the input data.  
To tackle this bias, we designed a method that evaluates the empirical significance of GO terms reported as enriched in modules. 

We used EMP to evaluate six AMI methods on GE and GWAS data. 
* jActiveModules (Ideker et al, 2002)
* NetBox (Cerami et al, 2010)
* HotNet2 (Raphael et al, 2015)
* KeyPathwayMiner (Alcaraz et al, 2012)
* Bionet (Beisser et al, 2010)
* [DOMINO](https://github.com/Shamir-Lab/DOMINO) - An algorithm we developed that produces reduced false GO term calls 

In our evaluation DOMINO outperformed the other algorithms.
  
A preprint version of the study is available at [BioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.10.984963v1)
  
If you wish to reproduce or extend our benchmark, please visit [the repository of the evaluation criteria](https://github.com/Shamir-Lab/EMP-benchmark).
## Outline


- [Set your environment](#set-your-environment)
- [Integrate your AMI algorithm with EMP](#integrate-your-ami-algorithm-with-emp)
- [Run EMP](#run-emp)
- [Main output files](#main-output-files)
- [EMP container](#emp-container)

## Set your environment

Download the sources and install according to the following:

Clone the repo from github:
```
git clone https://github.com/Shamir-Lab/EMP.git
cd EMP
```

EMP is written in Python 3.6. We recommend using a virtual environment. in Linux:
```
python -m venv emp-env
source emp-env/bin/activate
```

To install EMP dependencies type:
```
pip install -r  config/dependencies.txt
```


## Integrate your AMI algorithm with EMP

First, you need to make EMP to be aware to your AMI algorithm.    
1. Create an endpoint file your algorithm.  This file your algorithm should inherent AbstractRunner class (see `src/runnners/abstract_runners.py`).   
2. Make EMP acknowledged your algorithm by adding its instance to `ALGO_BY_NAMES` dictionary (under `src/runner/run_algo.py`). you can to it directly where `ALGO_BY_NAMES` is declared or at runtime   
 
For example see `src/emp/domino_runner.py` and `src/emp/run_algo.py`. 

## Run EMP

EMP consists of 6 main steps. For a specific set of input parameters, these steps should be carried sequentially.  
Each parameter can be specified as command line parameter (For example `python script.py --param1 value1 --param2 value2`). values of parameters which are not specified in the command line are taken from `config/conf.json`.      

1. `generate_permuted_datasets.py`: Generate permuted datasets.  
parameters:  
`--dataset_file`: path to dataset file.  
`--permuted_datasets_folder`: folder where permuted datasets reside.  
`--n_start`: starting positional index of permuted datasets.  
`--n_end`: ending positional index of permuted datasets.  
`--pf`: parallelization factor - number of cores EMP uses.  
`--override_permutations`: file of the biological network of the analysis.  
  
2. `generate_permuted_solutions.py`: Generate permuted solutions:  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--permuted_datasets_folder`: folder where permuted datasets reside.  
`--permuted_solutions_folder`: folder where permuted solutions reside.  
`--go_folder`: folder where GO files are located.  
`--n_start`: starting positional index of permuted datasets/solutions.  
`--n_end`: ending positional index of permuted datasets/solutions.  
`--pf`: parallelization factor - number of cores EMP uses.  
`--network_file`: file of the biological network of the analysis.  
`--override_permutations`: whether existing permutation with the same positional index should be deleted.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm.   
  
3. `generate_true_solution.py`: Generates an AMI solution based on the original (i.e. non-permuted) scores.  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--permuted_solutions_folder`: folder where permuted solutions reside.  
`--true_solutions_folder`: folder where true solutions reside.  
`--go_folder`: folder where GO files are located.  
`--network_file`: file of the biological network of the analysis.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm. 
  
4. `aggregare_bg_distribution.py`: Aggregates permuted solutions into a background distribution.  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--permuted_solutions_folder`: folder where permuted solutions reside.  
`--true_solutions_folder`: folder where true solutions reside.  
`--report_folder`: folder where analysis results reside.  
`--go_folder`: folder where GO files are located. This folder should contain the files "go.obo" GO term file, "gene2go" association file. The files for human are available at `http://purl.obolibrary.org/obo/go/go-basic.obo`, `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz`   
`--n_start`: starting positional index of permuted datasets/solutions.  
`--n_end`: ending positional index of permuted datasets/solutions.  
`--pf`: parallelization factor - number of cores EMP uses.  
  
5. `add_go_metadata.py`: Adds metadata for GO terms
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--go_folder`: folder where GO files are located. This folder should contain the files "go.obo" GO term file, "gene2go" association file. The files for human are available at `http://purl.obolibrary.org/obo/go/go-basic.obo`, `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz`   
`--report_folder`: folder where analysis results reside.  
`--n_permutations`: file of the biological network according which the analysis is carried.

6. `calculate_significance.py`: Calculates empirical p-values and correct for multiple testing (FDR):  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--go_folder`: folder where GO files are located. This folder should contain the files "go.obo" GO term file, "gene2go" association file. The files for human are available at `http://purl.obolibrary.org/obo/go/go-basic.obo`, `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz`  
`--report_folder`: folder where analysis results reside.  
`--n_total_samples`: enrichment score scores set size. This set is used to build the empirical distribution.  
`--n_dist_samples`: number of enrichment scores in the empirical distribution.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm.
  
Alternatively, you can run several steps sequentially with `run_emp.py`:


`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--permuted_datasets_folder`: folder where permuted datasets reside.  
`--permuted_solutions_folder`: folder where permuted solutions reside.  
`--true_solutions_folder`: folder where true solutions reside.  
`--report_folder`: folder where analysis results reside.  
`--go_folder`: folder where GO files are located. This folder should contain the files "go.obo" GO term file, "gene2go" association file. The files for human are available at `http://purl.obolibrary.org/obo/go/go-basic.obo`, `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz`.  
`--n_start`: starting positional index of permuted datasets/solutions.  
`--n_end`: ending positional index of permuted datasets/solutions.  
`--pf`: parallelization factor - number of cores EMP uses.  
`--network_file`: file of the biological network of the analysis. The file format should be *.sif with a single pair in each row. For example, see `data/emp_test/networks/dip.sif`    
`--override_permutations`: whether existing permutation with same positional index should be deleted.  
`--n_permutations`: file of the biological network according which the analysis is carried.  
`--n_total_samples`: enrichment score scores set size. This set is used to build the empirical distribution.  
`--n_dist_samples`: number of enrichment scores in the empirical distribution.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm  
`--processes`: comman delimited list of processes that should be carried.  names of processes are the same as its file name, excluding the file extension (.py) e.g "generate_permuted_datasets,generate_permuted_solution"  

Parameters default values are defined at `config/conf.json`  

## Main output files

`{permuted_datasets_folder}/{dataset_name}_random_{positional_index}`: Random datasets' folder. This is the output of `generate_permuted_datasets.py` script.   
`{permuted_solutions_folder}/sol_{algorithm}_{dataset_name}_random_{positional_index}`: Random solutions' folder. This is the output of `generate_permuted_solutions.py` script.  
`{true_solution_folder}/{dataset_name}_{algorithm}`: True solutions' folder. This is the output of `true_solution_folder.py` script.  
`{report_folder}/emp_diff_modules_{dataset_name}_{algorithm}.tsv`: aggregated scores file per GO term. This is the output of `aggregare_bg_distribution.py` script.  
`{report_folder}/emp_diff_modules_{dataset_name}_{algorithm}_md.tsv`: matrix of metadata for each GO term. This is the output of `add_go_metadata.py` script.  
`{report_folder}/emp_diff_modules_{dataset_name}_{algorithm}_passed_oob.tsv`: matrix of GO metadata + empirical p-value for each GO term. This is the final results of the pipeline. This is the output of `calculate_significance.py` script.  

fields of `emp_diff_modules_...` matrices:
* `GO id`: GO term id 
* `GO name`: GO term name
* `dist_n_samples`: List of GO terms' enrichment scores calculated on permuted datasets
* `hg_pval`: list of enrichment scores of each module per GO term of the real solution
* `hg_pval_max`: Max GO terms' enrichment scores over the real solution's modules
* `n_genes`: number of genes under each GO term (or its descendants) 
* `passed_oob_permutation_test`: number of genes under each GO term (or its descendants)
* `emp_pval`: list of empirical p-values of each module per GO term of the real solution
* `emp_pval_max`: number of genes under each GO term (or its descendants)
* `is_emp_pval_max_significant`: boolean indicating whether the empirical p-value is significant (i.e. below 0.05 th) 


## EMP container
EMP is also available as ready-to-use tool in a container
The container was generated and tested using udocker.It can also be loaded using Docker.
To load the container using udocker, do the following steps:
* Install [udocker](https://github.com/indigo-dc/udocker)
* Download the container from [here](https://drive.google.com/file/d/1kJLINKQNBMIUs-npavlzo8f4eOZf2Lls/view?usp=sharing)
* Load the container by running `udocker import --tocontainer --clone --name=emp emp-ubuntu-18.tar`
* Go inside the container by running `udocker run emp`
* the EMP project resides under /sandbox/
* EMP can be executed as described above in this README.
  * Advanced users can also be run EMP using /sandbox/main.sh.
    


