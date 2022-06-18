
bigwig_path=results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw
bias_bigwig=results/tobias/ATAC_PE/GM12878/data/GM12878/sorted_merged_corrected.bw

blacklist_region=$PWD/reference/GRch38_unified_blacklist.bed.gz
chrom_sizes=$PWD/reference/chrom.sizes
reference_fasta=$PWD/reference/hg38.genome.fa
output_dir=/mnt/lab_data2/anusri/chrombpnet/results/testing/
fold=$PWD/splits/fold_0.json
#chrom_sizes=/mnt/data/annotati

seed=5678
# defaults
inputlen=2114
outputlen=1000
filters=512
n_dilation_layers=8
negative_sampling_ratio=0.00000001

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

output_dir_n="results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/"
echo $( timestamp ): "python $PWD/src/training/test_generator.py \\
       --genome=$reference_fasta \\
       --bigwig=$bigwig_path \\
       --bias_bigwig=$bias_bigwig \\
       --peaks=$output_dir/filtered.peaks1.bed \\
       --params=$output_dir/chrombpnet_model_params.tsv \\
       --output_prefix=$output_dir/chrombpnet \\
       --chr_fold_path=$fold \\
       --batch_size=64 \\
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \\
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss" | tee -a $logfile
python $PWD/src/training/test_generator.py \
       --genome=$reference_fasta \
       --bigwig=$bigwig_path \
       --bias_bigwig=$bias_bigwig \
       --peaks=$output_dir/filtered.peaks1.bed \
       --params=$output_dir/chrombpnet_model_params.tsv \
       --output_prefix=$output_dir/chrombpnet \
       --chr_fold_path=$fold \
       --batch_size=64 \
       --architecture_from_file=$PWD/src/training/models/bpnet_model.py \
       --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss \
       --seed $seed

