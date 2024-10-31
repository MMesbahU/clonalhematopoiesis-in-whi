#!/bin/bash

## 

##
source scripts/useuse

use R-4.1

## input
my_rscipt=${1} # 

outDir=${2}

wl_csv=${3}

dnmt_csv=${4}

manrev_csv=${5}

out_prefix=${6}

all_vars=${7}

batch_name=${8}
##
out_allvarfile=${out_prefix}.${batch_name}.tsv

## 
Rscript ${my_rscipt} ${outDir} ${wl_csv} ${dnmt_csv} ${manrev_csv} ${out_prefix} ${all_vars} ${batch_name} ${out_allvarfile}

##

