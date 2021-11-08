#CUDA_VISIBLE_DEVICES=3 python score_snp_with_interpretation.py --output_path results/atac/ --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/unplug/model.0.hdf5
#python pr_curve_metrics.py --output_path results/atac/
CUDA_VISIBLE_DEVICES=2 python score_snp_with_interpretation.py --output_path results/atac_uncorrected/ --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/





