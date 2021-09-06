#!/bin/bash
module load meme

input_dir=$3
modisco_dir=$4
tomtom_dir=$5
meme_db=/oak/stanford/groups/akundaje/soumyak/motifs/motifs.meme.txt
dataset=$1


for mode in $2
do
	echo $mode
    [[ -d $tomtom_dir/$dataset ]] || mkdir $tomtom_dir/$dataset
    [[ -d $tomtom_dir/$dataset/$mode ]] || mkdir $tomtom_dir/$dataset/$mode

    for infile in $input_dir/$dataset/$mode/*
    do
        infile=$(basename $infile)

        for score_type in count_shap profile_shap
        do
            infile_path=$modisco_dir/$dataset/$mode/$infile.$score_type.hdf5
            if [[ -f $infile_path ]]
            then
                outfile_path=$tomtom_dir/$dataset/$mode/$infile.$score_type.tsv
                if [[ ! -f $outfile_path ]]
                then
                    echo $infile_path
		            python fetch_tomtom.py -m $infile_path -o $outfile_path -d $meme_db -n 10 &
                fi
            fi
        done
	done
done



