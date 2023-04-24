chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
#chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=$2
bias=$3
cellline=K562
gpu=1
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/negatives_data/negatives_subsample.bed
#regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/negatives_data/negatives_subsample.bed
output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/background_interpret/
#output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/background_interpret/




chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
#ref_fasta=/mnt/data/male.hg19.fa

#file=$output_dir/merged.$cellline
file=$output_dir/$cellline

#echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
#	-g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1


#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $bias -o $file --profile_or_counts profile

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts profile
python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1


#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $chrombpnet_nb -o $file --profile_or_counts profile
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/interpret/interpret.py -g $ref_fasta -r $regions -m $bias -o $file --profile_or_counts profile

#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions_v2.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1
#python src/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py -h5 $file.profile_scores.h5 -r $file.interpreted_regions.bed -c $chrom_sizes -o $file.profile.bw -s $file.profile.stat -t 1

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $output_dir/$cellline"bias.bw"  --output_path $file1"_scale_on_30k.txt"

