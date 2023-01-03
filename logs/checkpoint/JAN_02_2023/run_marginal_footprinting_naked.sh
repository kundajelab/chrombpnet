output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_05.13.2022_withnakedbias/chrombpnet_model/footprints_bias
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=chrombpnet_wo_bias.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5
model_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_05.13.2022_withnakedbias/chrombpnet_model/
gpu=0

mkdir $outout_dir

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir


