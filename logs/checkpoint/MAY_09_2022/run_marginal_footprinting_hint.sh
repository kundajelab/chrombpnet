output_dir=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=hint_atac.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
#motifs=AP1,CTCF

gpu=0

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting.sh $output_dir $reference_fasta $fold $model $motifs
