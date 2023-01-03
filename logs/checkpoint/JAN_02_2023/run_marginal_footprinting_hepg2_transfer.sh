output_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/chrombpnet_model/all_motifs_footprints_bias/
model_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/chrombpnet_model/

reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=chrombpnet_wo_bias.h5
motifs=GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1

gpu=0

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir
