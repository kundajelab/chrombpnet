chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
cellline=K562
gpu=0
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/data/30K.subsample.overlap.bed 
output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/dinuc_interpret_1/
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

mkdir $output_dir
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/dinuc_interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts counts

