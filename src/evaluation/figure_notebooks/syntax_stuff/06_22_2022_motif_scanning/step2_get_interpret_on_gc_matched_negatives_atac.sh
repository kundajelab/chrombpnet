gpu=3

# HEPG2 DNASE
out_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/interpret/
in_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/


ref_fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
negatives=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/negatives_data/negatives_with_summit.bed

shuf -n 5000 $negatives > $in_dir/5K.negatives.bed
 
regions=$in_dir/5K.negatives.bed

chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/nautilus_runs_jun16/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5 
file=$out_dir/5K.negatives

mode=profile
CUDA_VISIBLE_DEVICES=$gpu python ../../../../../src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts $mode

mode=counts
CUDA_VISIBLE_DEVICES=$gpu python ../../../../../src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts $mode


