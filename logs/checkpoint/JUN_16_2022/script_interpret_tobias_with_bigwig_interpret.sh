
python $PWD/src/evaluation/interpret/interpret.py \
       --genome=$PWD/reference/hg38.genome.fa \
       --regions=results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed \
       --output_prefix=$PWD/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_corrected/tobias_model/interpret/GM12878  \
       --model_h5=$PWD/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_corrected/tobias_model/tobias.h5

