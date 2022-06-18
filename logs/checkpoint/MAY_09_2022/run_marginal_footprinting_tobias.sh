output_dir=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=tobias_wo_bias.h5
motifs=GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
#motifs=AP1,CTCF

gpu=1

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting.sh $output_dir $reference_fasta $fold $model $motifs
