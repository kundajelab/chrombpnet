reference_fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
bigwig_path=$1
nonpeaks=$2
fold=$3
output_dir=$4
#logfile=$5

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

mkdir $output_dir/new_test_metrics/

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_dir"/new_test_metrics/test_metrics.log"
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
        --output_prefix=$output_dir/new_test_metrics/chrombpnet \\
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
        --output_prefix=$output_dir/new_test_metrics/chrombpnet \
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
        --output_prefix=$output_dir/new_test_metrics/chrombpnet_wo_bias \\
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
        --output_prefix=$output_dir/new_test_metrics/chrombpnet_wo_bias \
        --batch_size=256 \
        --model_h5=$output_dir/chrombpnet_wo_bias.h5 | tee -a $logfile

