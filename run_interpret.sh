biasdir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/H1ESC/H1ESC_07.17.2022_bias_128_4_1234_0.8_fold_4_data_type_ATAC_PE/bias_model/

ref_fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa

cell_line=H1ESC
CUDA_VISIBLE_DEVICES=1 python $PWD/src/evaluation/interpret/interpret.py \
        --genome=$ref_fasta \
        --regions=$biasdir/H1ESC.interpreted_regions.bed \
        --output_prefix=$biasdir/interpret/$cell_line \
        --model_h5=$biasdir/bias.h5
