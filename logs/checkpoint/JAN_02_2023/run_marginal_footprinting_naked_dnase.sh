output_dir=results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_05.13.2022_withnakedbias/chrombpnet_model/footprints_motifs_corrected/


mkdir $outout_dir
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=chrombpnet_wo_bias.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
model_dir=results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_05.13.2022_withnakedbias/chrombpnet_model/
gpu=1

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir
