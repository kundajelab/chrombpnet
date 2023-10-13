reference_fasta=$1
bigwig_path=$2
nonpeaks=$3
fold=$4
output_dir=$5
logfile=$6

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
    logfile=$output_dir"/test/test_metrics.log"
    touch $logfile
fi

inputlen=2114
outputlen=1000

# # # predictions and metrics on the chrombpnet model trained
echo $( timestamp ): "python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_new.py \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$output_dir/filtered.peaks.bed \\
        --nonpeaks=$nonpeaks \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/test/chrombpnet \\
        --batch_size=256 \\
        --model_h5=$output_dir/chrombpnet.h5" | tee -a $logfile
python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_new.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$nonpeaks \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/test/chrombpnet \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet.h5 | tee -a $logfile

# # # predictions and metrics on the chrombpnet model without bias trained
echo $( timestamp ): "python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_new.py \\
        --genome=$reference_fasta \\
        --bigwig=$bigwig_path \\
        --peaks=$output_dir/filtered.peaks.bed \\
        --nonpeaks=$nonpeaks \\
        --chr_fold_path=$fold \\
        --inputlen=$inputlen \\
        --outputlen=$outputlen \\
        --output_prefix=$output_dir/test/chrombpnet_wo_bias \\
        --batch_size=256 \\
        --model_h5=$output_dir/chrombpnet_wo_bias.h5 " | tee -a $logfile
python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_new.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$output_dir/filtered.peaks.bed \
        --nonpeaks=$nonpeaks \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/test/chrombpnet_wo_bias \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet_wo_bias.h5 | tee -a $logfile

