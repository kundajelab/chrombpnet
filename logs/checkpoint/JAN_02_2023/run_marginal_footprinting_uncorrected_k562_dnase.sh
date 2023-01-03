output_dir=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/footprints_motifs_uncorrected/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=hint_atac.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
#motifs=AP1,CTCF
model_dir=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/
mkdir $output_dir

gpu=7

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir

