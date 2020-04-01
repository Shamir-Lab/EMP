# EMP

EMP: EMpirical Pipeline for correcting GO terms obtained by network-based module discovery algorithms.

- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Main output files](#main-output-files)
- [Advanced usage](#advanced-usage)

## Installation

### From sources
Download the sources and install according to the following:

EMP is written in Python3. The necessary libraries will all be installed by the `setup.py` script.
We recommend using a virtual environment. For example, in Linux, before running `setup.py`:
```
python -m venv emp-env
source emp-env/bin/activate
```
To install, download and run setup.py:
```
    git clone https://github.com/hag007/EMP.git
    cd EMP
    python setup.py install
```
It is possible to install as a user without root permissions:
```
python setup.py install --user
```

## Basic Usage
You can run each step of EMP individually:

1. Generate permuted solutions:
```
generate_permuted_solutions --dataset </path/to/dataset> --algo <alg> --omic_type <ot> --out_dir</path/to/output_folder> \
--n_start <0> --n_end <5000> --pf <10>  # --override_permutations <False>
```

2. aggregate permuted solutions into background distribution, and create the real solution:
```
aggregate_bg_distribution --dataset </path/to/dataset> --algo <alg> --omic_type <ot> --out_dir</path/to/output_folder> \
--n_start <0> --n_end <5000> --pf <10> --calc_true_scores <true>
```

3. Add metadata for GO terms
```
add_go_metadata --dataset </path/to/dataset> --dataset </path/to/dataset> --algo <alg> --omic_type <ot> --out_dir</path/to/output_folder>
```

4. Calculate p-values and correct for multiple testing (FDR):
```
calculate_significance --dataset </path/to/dataset> --dataset </path/to/dataset> --algo <alg> --omic_type <ot> --out_dir</path/to/output_folder>
```

Alternatively, you can run the first step and then run step 2,3,4 together:
```
pipeline_bg_dist_to_pval --dataset </path/to/dataset> --algo <alg> --omic_type <ot> --out_dir</path/to/output_folder> \
--n_start <0> --n_end <5000> --pf <10> --calc_true_scores <true>
```


In order to run a particular algorithm, you need to create a "runner" file and add it to "runners" folder. This file should implement "runners/i_runner.py" class.
For example, see runners/domino_runner.py



## Main output files

