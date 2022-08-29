output_dir=$1
reference_fasta=$2
fold=$3
model=$4
motifs=$5
model_dir=$6

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

mkdir $output_dir

# create the log file
if [ -z "$logfile" ]
  then
    echo "No logfile supplied - creating one"
    logfile=$output_dir"/footprinting.log"
    touch $logfile
fi

output_file=$output_dir/motif

echo $( timestamp ): "mkdir $output_dir" | tee -a $logfile
        
echo $( timestamp ): "python $PWD/src/evaluation/marginal_footprints/marginal_footprinting_all.py \\
	-g $reference_fasta \\
	-r $model_dir/filtered.nonpeaks.bed \\
	--chr_fold_path=$fold \\
	-m $model_dir/$model \\
	-bs 256 \\
	-o $output_file \\
	-pwm_f src/evaluation/marginal_footprints/gm_benchmarking_motifs.tsv \\
	-mo $motifs" | tee -a $logfile

python $PWD/src/evaluation/marginal_footprints/marginal_footprinting_all.py \
	-g $reference_fasta \
	-r $model_dir/filtered.nonpeaks.bed \
	--chr_fold_path=$fold \
	-m $model_dir/$model \
	-bs 256 \
	-o $output_file \
	-pwm_f src/evaluation/marginal_footprints/gm_benchmarking_motifs.tsv \
	-mo $motifs | tee -a $logfile
