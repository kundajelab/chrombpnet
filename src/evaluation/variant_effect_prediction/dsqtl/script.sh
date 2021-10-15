CUDA_VISIBLE_DEVICES=3 python score_snp.py --output_path atac_output_uncorrected --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/model.0.hdf5
python pr_curve.py --output_path atac_output_uncorrected
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --output_path dnase_output_uncorrected --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/DNASE/GM12878/4_1_shifted_DNASE_10.04.2021_bias_filters_128/final_model_step3/model.0.hdf5
#python pr_curve.py --output_path dnase_output_uncorrected
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --output_path dnase_output1 --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_sept_09/GM12878/DNASE_08.11.2021/final_model_step3_new/unplug/model.0.hdf5
#python pr_curve.py --output_path dnase_output1
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --output_path dnase_output --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/DNASE/GM12878/4_1_shifted_DNASE_10.04.2021_bias_filters_128/final_model_step3/unplug/model.0.hdf5
#python pr_curve.py --output_path dnase_output
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --output_path atac_output --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/unplug/model.0.hdf5
#python pr_curve.py --output_path atac_output


