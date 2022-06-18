reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=chrombpnet.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1,GABPA,BACH1+MAFK,HNF4G,NFYB
output_dir=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/motifs_footprints_uncorrected/
model_dir=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/
mkdir $output_dir
gpu=0

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting.sh $output_dir $reference_fasta $fold $model $motifs $model_dir
