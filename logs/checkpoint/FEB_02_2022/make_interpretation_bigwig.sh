#regions=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/testing/beta_globin.coords
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/chr20.peaks.bed
#regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/chr20.peaks.bed
#regions=/srv/scratch/anusri/colloboration_data_overlap/interesting.locus


output_dir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/chr20_interpret_new/
moutput_dir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/
mkdir $output_dir
gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py \
        -bm $moutput_dir/bias_model_scaled.h5 \
        -cm $moutput_dir/chrombpnet.h5 \
        -cmb $moutput_dir/chrombpnet_wo_bias.h5 \
        -r $regions \
        -g /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
        -o $output_dir/pred \
        -c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -t 1 \



CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=$regions \
       --output_prefix=$output_dir/chrombpnet_wo_bias \
       --model_h5=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/chrombpnet_wo_bias.h5 \
       --profile_or_counts profile
#       --model_h5=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/bias_model_scaled.h5 \



python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
	-h5 $output_dir/chrombpnet_wo_bias.profile_scores.h5 \
	-o $output_dir/chrombpnet_wo_bias_profile.bw \
	-s $output_dir/chrombpnet_wo_bias_profile.stats \
	-t 1


CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=$regions \
       --output_prefix=$output_dir/chrombpnet_wo_bias \
       --model_h5=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/chrombpnet_wo_bias.h5 \
       --profile_or_counts counts
#       --model_h5=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/bias_model_scaled.h5 \

python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
        -r $regions \
	-h5 $output_dir/chrombpnet_wo_bias.counts_scores.h5 \
	-o $output_dir/chrombpnet_wo_bias_counts.bw \
	-s $output_dir/chrombpnet_wo_bias_counts.stats \
	-t 1



#CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
#       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
#       --regions=$regions \
#       --output_prefix=$output_dir/chrombpnet \
#       --model_h5=results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022_bigger_jitter_on_bias/chrombpnet_model/chrombpnet.h5 \



#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
#	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
#        -r $regions \
#	-h5 $output_dir/chrombpnet.profile_scores.h5 \
#	-o $output_dir/chrombpnet_profile.bw \
#	-s $output_dir/chrombpnet_profile.stats \
#	-t 1

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py \
#	-c /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
#        -r $regions \
#	-h5 $output_dir/chrombpnet.counts_scores.h5 \
#	-o $output_dir/chrombpnet_counts.bw \
#	-s $output_dir/chrombpnet_counts.stats \
#	-t 1


