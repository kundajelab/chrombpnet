#!/bin/bash

in_bam=$1
bigwig_prefix=$2
data_type=$3
reference_fasta=$4
chrom_sizes=$5

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

logfile=$bigwig_prefix"_preprocessing.log"
touch $logfile

bash $PWD/src/helpers/preprocessing/bam_to_bigwig.sh $in_bam $bigwig_prefix $data_type $chrom_sizes $logfile

echo $( timestamp ): "python $PWD/src/helpers/preprocessing/analysis/build_pwm_from_bigwig.py -i $bigwig_prefix_unstranded.bw -g $reference_fasta -o $bigwig_prefix_bias_pwm.png -c chr20 -cz $chrom_sizes" | tee -a $logfile
python $PWD/src/helpers/preprocessing/analysis/build_pwm_from_bigwig.py -i $bigwig_prefix"_unstranded.bw" -g $reference_fasta -o $bigwig_prefix"_bias_pwm.png" -c "chr20" -cz $chrom_sizes 
