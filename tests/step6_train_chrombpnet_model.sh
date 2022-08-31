#!/bin/bash

echo "WARNING: If upgrading from v1.0 or v1.1 to v1.2. Note that chrombpnet has undergone linting to generate a modular structure for release on pypi.Hard-coded script paths are no longer necessary. Please refer to the updated README (below) to ensure your script calls are compatible with v1.2"

# exit when any command fails
set -e
set -o pipefail

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

reference_fasta=${1?param missing - reference_fasta}
bigwig_path=${2?param missing - bigwig_path }
overlap_peak=${3?param missing - overlap_peak}
nonpeaks=${4?param missing - nonpeaks}
fold=${5?param missing - fold}
bias_model=${6?param missing - bias_model}
output_dir=${7?param missing - output_dir}
data_type=${8?param missing - data_type}
seed=${9:-1234}
logfile=${10} #optional
pwm_f=${11}  #optional

# defaults
inputlen=2114
outputlen=1000
filters=512
n_dilation_layers=8
negative_sampling_ratio=0.1

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_dir"/train_chrombpnet_model.log"
    touch $logfile
fi



# marginal footprinting
if [[ "$data_type" = "DNASE_SE" || "$data_type" = "DNASE_PE" ]] ; then
    echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
    mkdir $output_dir/footprints
            echo $( timestamp ): "chrombpnet_marginal_footprints \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        -m $output_dir/chrombpnet_wo_bias.h5 \\
        -bs 512 \\
        -o $output_dir/footprints/corrected \\
        -pwm_f $pwm_f "| tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/chrombpnet_wo_bias.h5 \
		-bs 512 \
		-o $output_dir/footprints/corrected \
		-pwm_f $pwm_f | tee -a $logfile
elif [[ "$data_type" = "ATAC_SE" || "$data_type" = "ATAC_PE"  ]] ; then
    echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
    mkdir $output_dir/footprints
            echo $( timestamp ): "chrombpnet_marginal_footprints \\
        -g $reference_fasta \\                     
        -r $output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        -m $output_dir/chrombpnet_wo_bias.h5 \\
        -bs 512 \\
        -o $output_dir/footprints/corrected \\
        -pwm_f $pwm_f  | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/chrombpnet_wo_bias.h5 \
		-bs 512 \
		-o $output_dir/footprints/corrected \
		-pwm_f $pwm_f | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi

# marginal footprinting bias model
if [[ "$data_type" = "DNASE_SE" || "$data_type" = "DNASE_PE" ]] ; then
    echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
    mkdir $output_dir/footprints
            echo $( timestamp ): "chrombpnet_marginal_footprints \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        -m $output_dir/bias_model_scaled.h5 \\
        -bs 512 \\
        -o $output_dir/footprints/bias \\
        -pwm_f $pwm_f  | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/bias_model_scaled.h5 \
		-bs 512 \
		-o $output_dir/footprints/bias \
		-pwm_f $pwm_f | tee -a $logfile
elif [[ "$data_type" = "ATAC_SE" || "$data_type" = "ATAC_PE"  ]] ; then
    echo $( timestamp ): "mkdir $output_dir/footprints" | tee -a $logfile
    mkdir $output_dir/footprints
    echo $( timestamp ): "chrombpnet_marginal_footprints \\
    	 -g $reference_fasta \\
         -r $output_dir/filtered.nonpeaks.bed \\
         -chr "chr1" \\
         -m $output_dir/bias_model_scaled.h5 \\
         -bs 512 \\
         -o $output_dir/footprints/bias \\
         -pwm_f $pwm_f" | tee -a $logfile
    chrombpnet_marginal_footprints \
	-g $reference_fasta \
	-r $output_dir/filtered.nonpeaks.bed \
	--chr_fold_path=$fold \
	-m $output_dir/bias_model_scaled.h5 \
	-bs 512 \
	-o $output_dir/footprints/bias \
	-pwm_f $pwm_f  | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi

