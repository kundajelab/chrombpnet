#!/bin/bash

ml python/3.6.1
ml py-numpy/1.18.1_py36
export PYTHONPATH=$PYTHONPATH:/home/users/anusri/.local/lib/python3.6/site-packages/

input_dir=$3
output_dir=$4
dataset=$1
cores=20
seqlets=40000

for mode in $2
do
    echo $input_dir/$dataset/$mode/
    [[ -d $output_dir/$dataset ]] || mkdir $output_dir/$dataset
    [[ -d $output_dir/$dataset/$mode ]] || mkdir $output_dir/$dataset/$mode

    for infile in $input_dir/$dataset/$mode/*
    do
        infile=$(basename $infile)
        infile_path=$input_dir/$dataset/$mode/$infile

        for score_type in count_shap profile_shap
        do
            outfile_path=$output_dir/$dataset/$mode/$infile.$score_type

            if [[ ! -f $outfile_path.hdf5 ]]
            then
                echo $infile_path
                echo $outfile_path

                sbatch --export=ALL --requeue \
                       -J $mode.$infile \
                       -p owners,akundaje \
                       -t 720 -c $cores --mem=150G \
                       -o $outfile_path.log.o \
                       -e $outfile_path.log.e \
                       run_modisco.sh $infile_path $outfile_path.hdf5 $score_type $cores $seqlets
            fi
        done
    done
done

