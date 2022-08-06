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


# this script does the following -  
# (1) filters your peaks/nonpeaks (removes outliers and removes edge cases and creates a new filtered set)
# (2) scales the given bias model on the non-peaks
# (3) Calculates the counts loss weight 
# (4) Creates a TSV file that can be loaded into the next step
echo $( timestamp ): "chrombpnet_hyperparams \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$overlap_peak \\
       --nonpeaks=$nonpeaks \\
       --outlier_threshold=0.99 \\
       --chr_fold_path=$fold \\
       --negative_sampling_ratio=$negative_sampling_ratio \\
       --inputlen=$inputlen \\
       --outputlen=$outputlen \\
       --max_jitter=500 \\
       --filters=$filters \\
       --n_dilation_layers=$n_dilation_layers \\
       --bias_model_path=$bias_model \\
       --output_dir=$output_dir " | tee -a $logfile
chrombpnet_hyperparams \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$overlap_peak \
    --nonpeaks=$nonpeaks \
    --outlier_threshold=0.99 \
    --chr_fold_path=$fold \
    --negative_sampling_ratio=$negative_sampling_ratio \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --max_jitter=500 \
    --filters=$filters \
    --n_dilation_layers=$n_dilation_layers \
    --bias_model_path=$bias_model \
    --output_dir=$output_dir | tee -a $logfile

## this script does the following -  
# (1) trains a model on the given peaks/nonpeaks
# (2) The parameters file input to this script should be TSV seperatedp
chrombpnet_with_bias_model_path=`which chrombpnet_with_bias_model.py`
echo $( timestamp ): "chrombpnet_train \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --peaks=$output_dir/filtered.peaks.bed \\
       --nonpeaks=$output_dir/filtered.nonpeaks.bed \\
       --params=$output_dir/chrombpnet_model_params.tsv \\
       --output_prefix=$output_dir/chrombpnet \\
       --chr_fold_path=$fold \\
       --seed=$seed \\
       --batch_size=64 \\
       --architecture_from_file=$chrombpnet_with_bias_model_path \\
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss" | tee -a $logfile

chrombpnet_train \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$output_dir/filtered.peaks.bed \
    --nonpeaks=$output_dir/filtered.nonpeaks.bed \
    --params=$output_dir/chrombpnet_model_params.tsv \
    --output_prefix=$output_dir/chrombpnet \
    --chr_fold_path=$fold \
    --seed=$seed \
    --batch_size=64 \
    --architecture_from_file=$chrombpnet_with_bias_model_path \
    --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss | tee -a $logfile

#predictions and metrics on the chrombpnet model trained
echo $( timestamp ): "chrombpnet_predict \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$output_dir/filtered.peaks.bed \\
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/chrombpnet \\
        --batch_size=256 \\
        --model_h5=$output_dir/chrombpnet.h5" | tee -a $logfile
chrombpnet_predict \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$output_dir/filtered.peaks.bed \
    --nonpeaks=$output_dir/filtered.nonpeaks.bed \
    --chr_fold_path=$fold \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --output_prefix=$output_dir/chrombpnet \
    --batch_size=256 \
    --model_h5=$output_dir/chrombpnet.h5 | tee -a $logfile

# predictions and metrics on the chrombpnet model without bias trained
echo $( timestamp ): "chrombpnet_predict \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$output_dir/filtered.peaks.bed \\
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/chrombpnet_wo_bias \\
        --batch_size=256 \\
        --model_h5=$output_dir/chrombpnet_wo_bias.h5 " | tee -a $logfile
chrombpnet_predict \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$output_dir/filtered.peaks.bed \
    --nonpeaks=$output_dir/filtered.nonpeaks.bed \
    --chr_fold_path=$fold \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --output_prefix=$output_dir/chrombpnet_wo_bias \
    --batch_size=256 \
    --model_h5=$output_dir/chrombpnet_wo_bias.h5 | tee -a $logfile

#predictions and metrics on the bias model trained
echo $( timestamp ): "chrombpnet_predict \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$output_dir/filtered.peaks.bed \\
        --nonpeaks=$output_dir/filtered.nonpeaks.bed \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/bias \\
        --batch_size=256 \\
        --model_h5=$output_dir/bias_model_scaled.h5" | tee -a $logfile
chrombpnet_predict \
    --genome=$reference_fasta \
    --bigwig=$bigwig_path \
    --peaks=$output_dir/filtered.peaks.bed \
    --nonpeaks=$output_dir/filtered.nonpeaks.bed \
    --chr_fold_path=$fold \
    --inputlen=$inputlen \
    --outputlen=$outputlen \
    --output_prefix=$output_dir/bias \
    --batch_size=256 \
    --model_h5=$output_dir/bias_model_scaled.h5 | tee -a $logfile

# marginal footprinting



#path to pwm file
if [ -z "$pwm_f" ]
  then
      echo "No pwm file supplied, using default"
      tee -a motif_to_pwm.default.tsv <<EOF
EOF 
tn5_1	GCACAGTACAGAGCTG
tn5_2	GTGCACAGTTCTAGAGTGTGCAG
tn5_3	CCTCTACACTGTGCAGAA
tn5_4	GCACAGTTCTAGACTGTGCAG
tn5_5	CTGCACAGTGTAGAGTTGTGC
dnase_1	TTTACAAGTCCA
dnase_2	TGTACTTACGAA
NRF1	GCGCATGCGC
AP1	CGATATGACTCATCCC
CTCF	TTGGCCACTAGGGGGCGCTAT
ETS	CCGAAAGCGGAAGTGAGAC	
SP1	AAGGGGGCGGGGCCTAA
RUNX	CCCTAACCACAGCCC
NFKB	GCAAGGGAAATTCCCCAGG
GATA+TAL	GGCTGGGGGGGGCAGATAAGGCC
TAL	GGCTGGG
NFYB	CCAGCCAATCAGAGC
GABPA	GAAACCGGAAGTGGCC
BACH1+MAFK	AACTGCTGAGTCATCCCG
NRF1	CCCCGCGCATGCGCAGTGC
HNF4G	CCGTTGGACTTTGGACCCTG
EOF
      pwm_f=motif_to_pwm.default.tsv
fi

EOF

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
        -pwm_f $pwm_f \\
        -mo dnase_1,dnase_2" | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/chrombpnet_wo_bias.h5 \
		-bs 512 \
		-o $output_dir/footprints/corrected \
		-pwm_f $pwm_f \
		-mo dnase_1,dnase_2 | tee -a $logfile
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
        -pwm_f $pwm_f \\
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5" | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/chrombpnet_wo_bias.h5 \
		-bs 512 \
		-o $output_dir/footprints/corrected \
		-pwm_f $pwm_f \
		-mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5 | tee -a $logfile
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
        -pwm_f $pwm_f \\
        -mo dnase_1,dnase_2" | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/bias_model_scaled.h5 \
		-bs 512 \
		-o $output_dir/footprints/bias \
		-pwm_f $pwm_f \
		-mo dnase_1,dnase_2 | tee -a $logfile
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
        -pwm_f $pwm_f \\
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5" | tee -a $logfile
	    chrombpnet_marginal_footprints \
		-g $reference_fasta \
		-r $output_dir/filtered.nonpeaks.bed \
		--chr_fold_path=$fold \
		-m $output_dir/bias_model_scaled.h5 \
		-bs 512 \
		-o $output_dir/footprints/bias \
		-pwm_f $pwm_f \
		-mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5 | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi

