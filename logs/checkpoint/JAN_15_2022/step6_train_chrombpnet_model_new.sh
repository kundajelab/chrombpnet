reference_fasta=$1
bigwig_path=$2
overlap_peak=$3
nonpeaks=$4
fold=$5
bias_model=$6
output_dir=$7
data_type=$8
logfile=$9

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

# marginal footprinting - for corrected model

if [[ "$data_type" = "DNASE_SE" || "$data_type" = "DNASE_PE" ]] ; then
        echo $( timestamp ): "mkdir $output_dir/dnase_footprints" | tee -a $logfile
        mkdir $output_dir/dnase_footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        -chr "chr1" \\
        -m $output_dir/chrombpnet_wo_bias.h5 \\
        -bs 512 \\
        -o $output_dir/dnase_footprints/corrected \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo dnase_1,dnase_2" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet_wo_bias.h5 \
        -bs 512 \
        -o $output_dir/dnase_footprints/corrected \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo dnase_1,dnase_2 | tee -a $logfile
elif [[ "$data_type" = "ATAC_SE" || "$data_type" = "ATAC_PE"  ]] ; then
        echo $( timestamp ): "mkdir $output_dir/tn5_footprints" | tee -a $logfile
        mkdir $output_dir/tn5_footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        -chr "chr1" \\
        -m $output_dir/chrombpnet_wo_bias.h5 \\
        -bs 512 \\
        -o $output_dir/tn5_footprints/corrected \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet_wo_bias.h5 \
        -bs 512 \
        -o $output_dir/tn5_footprints/corrected \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5 | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi


# marginal footprtining bias model
if [[ "$data_type" = "DNASE_SE" || "$data_type" = "DNASE_PE" ]] ; then
        echo $( timestamp ): "mkdir $output_dir/dnase_footprints" | tee -a $logfile
        mkdir $output_dir/dnase_footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        -chr "chr1" \\
        -m $output_dir/bias_model_scaled.h5 \\
        -bs 512 \\
        -o $output_dir/dnase_footprints/bias \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo dnase_1,dnase_2" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/bias_model_scaled.h5 \
        -bs 512 \
        -o $output_dir/dnase_footprints/bias \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo dnase_1,dnase_2 | tee -a $logfile
elif [[ "$data_type" = "ATAC_SE" || "$data_type" = "ATAC_PE"  ]] ; then
        echo $( timestamp ): "mkdir $output_dir/tn5_footprints" | tee -a $logfile
        mkdir $output_dir/tn5_footprints
        echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \\
        -g $reference_fasta \\
        -r $output_dir/filtered.nonpeaks.bed \\
        -chr "chr1" \\
        -m $output_dir/bias_model_scaled.h5 \\
        -bs 512 \\
        -o $output_dir/tn5_footprints/bias \\
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \\
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5" | tee -a $logfile
        python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/bias_model_scaled.h5 \
        -bs 512 \
        -o $output_dir/tn5_footprints/bias \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5 | tee -a $logfile
else
    echo "ERROR: unknown data type " $data_type | tee -a $logfile
fi

