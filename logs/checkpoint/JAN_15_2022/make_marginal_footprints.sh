cell_line=K562
data_type="DNASE_SE"
setting=DNASE_SE_12.30.2021
gpu=3

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
        -mo NFKB,GATA+TAL,CTCF,GABPA,BACH1+MAFK,NFYB,HNF4G
#        -mo NRF1,CTCF,ETS,SP1,RUNX,NFKB,GATA+TAL,dnase_1,dnase_2

CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/marginal_footprints/marginal_footprinting.py \
        -g $reference_fasta \
        -r $output_dir/filtered.nonpeaks.bed \
        -chr "chr1" \
        -m $output_dir/chrombpnet.h5 \
        -bs 256 \
        -o $output_dir/footprints/chrombpnet \
        -pwm_f src/evaluation/marginal_footprints/motif_to_pwm.tsv \
        -mo NFKB,GATA+TAL,CTCF,GABPA,BACH1+MAFK,NFYB,HNF4G
#        -mo NRF1,CTCF,ETS,SP1,RUNX,NFKB,GATA+TAL,dnase_1,dnase_2
