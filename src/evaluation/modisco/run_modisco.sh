#!/bin/bash

ml python/3.6.1
ml py-numpy/1.18.1_py36
export PYTHONPATH=$PYTHONPATH:/home/users/anusri/.local/lib/python3.6/site-packages/

infile=$1
outfile=$2
score_type=$3
cores=$4
seqlets=$5

python3 run_modisco.py $infile $outfile $score_type $cores $seqlets

