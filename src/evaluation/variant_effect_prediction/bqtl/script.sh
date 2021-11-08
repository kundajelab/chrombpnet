#CUDA_VISIBLE_DEVICES=1 python score_snp.py --subsample --output_path testing --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/DNASE/GM12878/4_1_shifted_DNASE_10.04.2021_bias_filters_128/final_model_step3/unplug/model.0.hdf5
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --subsample --output_path dnase_subsample --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/DNASE/GM12878/4_1_shifted_DNASE_10.04.2021_bias_filters_128/final_model_step3/unplug/model.0.hdf5
#CUDA_VISIBLE_DEVICES=3 python score_snp.py --subsample --output_path atac_subsample --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/unplug/model.0.hdf5
#CUDA_VISIBLE_DEVICES=3 python score_snp_with_interpretation.py --subsample --output_path results/atac_uncorrected --model_file /srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/
CUDA_VISIBLE_DEVICES=3 python score_snp_with_interpretation.py --subsample --output_path results/spi1_chip --model_file /srv/scratch/annashch/chrombpnet/gm12878_tf_dnase_histone/Spi1ProfileModel_RevComp_9layers.h5



