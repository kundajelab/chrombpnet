chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet.h5
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/bias_model_scaled.h5
cellline=K562
gpu=0


fold=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json
observed=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/K562_unstranded.bw
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/data/peaks_no_blacklist.bed
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/ATAC/K562/scale_fg/isotonic/ENCSR868FGK/
mkdir $output_dir


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_normalized_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_normalized_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -obs $observed -f $fold -t 1


