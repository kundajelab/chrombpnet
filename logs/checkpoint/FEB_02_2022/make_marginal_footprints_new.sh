cell_line=GM12878
data_type="DNASE_SE"
#setting=ATAC_PE_12.30.2021
#setting=DNASE_SE_withk562bias_01.03.2022
setting=DNASE_SE_01.15.2022_0.9
gpu=0

reference_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
output_dir=$main_dir/$setting/chrombpnet_model/
CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet_wo_bias.h5 \
        -bs 256 \
        -o $output_dir/footprints/chrombpnet_wo_bias \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo NRF1,CTCF,ETS,SP1,RUNX,NFKB,GATA+TAL,dnase_1,dnase_2,AP1

CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet.h5 \
        -bs 256 \
        -o $output_dir/footprints/chrombpnet \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo NRF1,CTCF,ETS,SP1,RUNX,NFKB,GATA+TAL,dnase_1,dnase_2,AP1
