
gpu=1
ref_fasta=$PWD/reference/hg38.genome.fa
regions=/mnt/lab_data2/anusri/chrombpnet/results/print/ATAC_PE/ENCSR158XTU/chrombpnet_compare/ENCSR158XTU/auxiliary/30K_subsample_peaks.bed
output_dir=/mnt/lab_data2/anusri/chrombpnet/results/print/ATAC_PE/ENCSR158XTU/ATAC_PE_05.09.2025_print_with_bias_bigwig/ 
cell_line=print_ENCSR158XTU

mkdir=$output_dir/print_model/interpret_new/

CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$regions \
        --output_prefix=$output_dir/print_model/interpret_new/$cell_line \
        --model_h5=$output_dir/print_model/print_wo_bias.h5  | tee -a $logfile

