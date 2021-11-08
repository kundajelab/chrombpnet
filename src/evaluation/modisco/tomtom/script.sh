flank=500
input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/SIGNAL/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/SIGNAL/

#model_name=4_4_shifted_ATAC_09.06.2021_bias_filters_500_withgm12878bias
#model_name=4_4_shifted_ATAC_09.16.2021_subsample_5M_gm12878biasfit
#model_name=4_4_shifted_ATAC_09.16.2021_subsample_100M_k562biasfit
#model_name=4_4_shifted_ATAC_09.21.2021_bias_filters_128_new
#model_name=4_4_shifted_ATAC_09.29.2021_bias_filters_500_new
#model_name=4_4_shifted_ATAC_09.30.2021_subsample_100M_k562biasfit
#model_name=4_4_shifted_ATAC_10.01.2021_subsample_50M_k562biasfit
#model_name=4_1_shifted_DNASE_10.04.2021_bias_filters_128_new
#model_name=4_4_shifted_ATAC_10.05.2021_bias_filters_500subsample_50M_new
#model_name=ATAC_10.09.2021_withuniversalbias_universalbiasfit
#model_name=ATAC_10.09.2021_withinvivobias_new
#model_name=ATAC_10.14.2021_withuniversalbias_500filt_universalbiasfit
#model_name=ATAC_10.14.2021_withinvivobias_new
model_name=ATAC_10.14.2021_withinvivobias_500filts_mincount_new

#bash run_tomtom.sh SURAG $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh GM12878 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh HEPG2 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh IMR90 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh H1 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh K562 $model_name $input_dir $output_dir $tomtom_dir


input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/BIAS/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/BIAS/

#model_name=4_4_shifted_ATAC_09.06.2021_bias_filters_500
#model_name=4_4_shifted_ATAC_09.21.2021_bias_filters_128
#model_name=4_4_shifted_ATAC_09.29.2021_bias_filters_500
#model_name=4_4_shifted_ATAC_09.30.2021_bias_filters_128
#model_name=4_4_shifted_ATAC_10.04.2021_bias_filters_128
#model_name=4_1_shifted_DNASE_10.04.2021_bias_filters_128
#model_name=ATAC_10.09.2021_withinvivobias
#model_name=4_1_shifted_DNASE_10.04.2021_bias_filters_128
#model_name=ATAC_10.14.2021_withinvivobias
model_name=ATAC_10.14.2021_withinvivobias_500filts_mincount

#bash run_tomtom.sh naked_dna $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh GM12878 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh HEPG2 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh IMR90 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh H1 $model_name $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh K562 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh UNIVERSAL $model_name $input_dir $output_dir $tomtom_dir


