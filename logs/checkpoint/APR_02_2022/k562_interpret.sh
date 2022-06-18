#regions=results/chrombpnet/ATAC_PE/K562/data/30K.subsample.overlap.bed
#regions=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/testing/beta_globin.coords
regions=/srv/scratch/anusri/colloboration_data_overlap/interesting.locus

mkdir results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_256_8_2356/bias_model/interpret_locus
output_dir=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_256_8_2356/bias_model/interpret_locus
CUDA_VISIBLE_DEVICES=3 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/interpret/interpret.py  \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=$regions  \
       --output_prefix=$output_dir/interesting.locus    \
       --model_h5=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_256_8_2356/bias_model/bias.h5 \

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
        -c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
        -h5 $output_dir/interesting.locus.profile_scores.h5 \
        -o $output_dir/interesting.locus_profile.bw \
        -s $output_dir/interesting.locus_profile.stats \
        -t 1

