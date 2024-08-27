cellline=$1
biasth=$2
foldn=$3
outputdir=$4
ddatype=$5
echo "python /mnt/lab_data2/anusri/chrombpnet/src/helpers/hyperparameters/find_bias_hyperparams.py \\
       --genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \\
       --bigwig=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/data/$cellline"_unstranded.bw" \\
       --peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/data/peaks_no_blacklist.bed \\
       --nonpeaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/negatives_data_$foldn/negatives_with_summit.be>
       --outlier_threshold=0.99 \\
       --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_$foldn.json \\
       --inputlen=2114 \\
       --outputlen=1000 \\
       --max_jitter=0 \\
       --filters=128 \\
       --n_dilation_layers=4 \\
       --bias_threshold_factor=$biasth \\
       --output_dir $outputdir"

python /mnt/lab_data2/anusri/chrombpnet/src/helpers/hyperparameters/find_bias_hyperparams.py \
       --genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
       --bigwig=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/data/$cellline"_unstranded.bw" \
       --peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/data/peaks_no_blacklist.bed \
       --nonpeaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$ddatype/$cellline/negatives_data_$foldn/negatives_with_summit.bed \
       --outlier_threshold=0.99 \
       --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_$foldn.json \
       --inputlen=2114 \
       --outputlen=1000 \
       --max_jitter=0 \
       --filters=128 \
       --n_dilation_layers=4 \
       --bias_threshold_factor=$biasth \
       --output_dir $outputdir

