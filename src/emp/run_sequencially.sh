#!/bin/bash

pf=15

source ../../python/bin/activate
declare -a algs=("jactivemodules_greedy")
# declare -a dss=("ers" "iem" "apo" "cbx" "ift")
declare -a dss=("bmd" "hgt") #  "hc" "ror" "shera" "shezh")
#declare -a dss=("brca" "crh" "scz" "tri" "t2d") 
# declare -a dss=("bmd" "hgt" "cad" "amd" "af") 
# declare -a dss=("tnfa" "hc" "ror" "shera" "shezh" "ers" "iem" "apo" "cbx" "ift")
# declare -a dss=("brca" "crh" "scz" "tri" "t2d" "bmd" "hgt" "cad" "amd" "af")

script=$1
algs=(${3//,/ })
dss=(${2//,/ })
pf=${4}
 
for alg in "${algs[@]}" 
    do 
    for ds in "${dss[@]}" 
        do
        python ${1}.py --algo ${alg} --dataset_file /specific/netapp5/gaga/hagailevi/emp_test/original_datasets/${ds}.tsv  --pf ${pf}
    done
done

