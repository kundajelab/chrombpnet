output_dir=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.08.2022_tobias_corrected_not_softmax_custom_shift/tobias_model/footprints_motifs/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=tobias.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
#motifs=AP1,CTCF
model_dir=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.08.2022_tobias_corrected_not_softmax_custom_shift/tobias_model/

gpu=0

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting.sh $output_dir $reference_fasta $fold $model $motifs $model_dir
