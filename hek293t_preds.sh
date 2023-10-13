

chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/DNASE_SE_03.16.2023_1.0/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/DNASE_SE_03.16.2023_1.0/chrombpnet_model/chrombpnet.h5
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/DNASE_SE_03.16.2023_1.0/chrombpnet_model/bias_model_scaled.h5
cellline=ENCSR000EJR
gpu=0


regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/ENCSR000EJR/data/peaks_no_blacklist.bed
output_dir=/oak/stanford/groups/akundaje/anusri/hek293t_preds/


chrom_sizes=$PWD/reference/chrom.sizes
ref_fasta=$PWD/reference/hg38.genome.fa
file=$output_dir/$cellline

echo "CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1"
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/make_bigwigs/predict_to_bigwig.py -bm $bias -cm $chrombpnet -cmb $chrombpnet_nb --regions $regions \
        -g $ref_fasta -c $chrom_sizes -o $output_dir/$cellline -t 1
