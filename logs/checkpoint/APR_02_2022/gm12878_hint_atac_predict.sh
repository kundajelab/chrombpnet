bigwig_path=results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw
overlap_peak=$PWD/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed
blacklist_region=$PWD/reference/GRch38_unified_blacklist.bed.gz
chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
fold=$PWD/splits/fold_0.json
inputlen=2114
outputlen=1000
output_dir=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/

CUDA_VISIBLE_python $PWD/src/training/predict.py \
        --genome=$reference_fasta \
        --bigwig=$bigwig_path \
        --peaks=$overlap_peak \
        --chr_fold_path=$fold \
        --inputlen=$inputlen \
        --outputlen=$outputlen \
        --output_prefix=$output_dir/hint_atac_orign \
        --batch_size=256 \
        --model_h5=$output_dir/hint_atac.h5 | tee -a $logfile

