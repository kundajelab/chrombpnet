fold_0=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0
fold_1=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.08.2022_bias_128_4_1234_0.4_fold_1_data_type_ATAC_PE
fold_2=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.08.2022_bias_128_4_1234_0.4_fold_2_data_type_ATAC_PE
fold_3=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.14.2022_bias_128_4_1234_0.4_fold_3_data_type_ATAC_PE
fold_4=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.07.2022_bias_128_4_1234_0.4_fold_4_data_type_ATAC_PE

split_0=splits/fold_0.json
split_1=splits/fold_1.json
split_2=splits/fold_2.json
split_3=splits/fold_3.json
split_4=splits/fold_4.json

fold_0_model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_250M/GM12878_250M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE
fold_1_model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_250M/GM12878_250M_07.18.2022_bias_transfer_1234_fold_1_data_type_ATAC_PE
fold_2_model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_250M/GM12878_250M_07.19.2022_bias_transfer_1234_fold_2_data_type_ATAC_PE
fold_3_model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_250M/GM12878_250M_07.18.2022_bias_transfer_1234_fold_3_data_type_ATAC_PE
fold_4_model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_250M/GM12878_250M_07.18.2022_bias_transfer_1234_fold_4_data_type_ATAC_PE


CUDA_VISIBLE_DEVICES=0 bash predict_wrapper.sh $fold_0_model_dir $split_0 $fold_0/chrombpnet_model 
bash predict_wrapper.sh $fold_0_model_dir $split_1 $fold_1/chrombpnet_model 
bash predict_wrapper.sh $fold_0_model_dir $split_2 $fold_2/chrombpnet_model 
bash predict_wrapper.sh $fold_0_model_dir $split_3 $fold_3/chrombpnet_model 
bash predict_wrapper.sh $fold_0_model_dir $split_4 $fold_4/chrombpnet_model 



fold_0,100M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE
fold_1,100M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.18.2022_bias_transfer_1234_fold_1_data_type_ATAC_PE
fold_2,100M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.18.2022_bias_transfer_1234_fold_2_data_type_ATAC_PE
fold_3,100M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.18.2022_bias_transfer_1234_fold_3_data_type_ATAC_PE
fold_4,100M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.18.2022_bias_transfer_1234_fold_4_data_type_ATAC_PE
fold_0,50M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_50M/GM12878_50M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE
fold_1,50M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_50M/GM12878_50M_07.18.2022_bias_transfer_1234_fold_1_data_type_ATAC_PE
fold_2,50M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_50M/GM12878_50M_07.18.2022_bias_transfer_1234_fold_2_data_type_ATAC_PE
fold_3,50M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_50M/GM12878_50M_07.18.2022_bias_transfer_1234_fold_3_data_type_ATAC_PE
fold_4,50M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_50M/GM12878_50M_07.19.2022_bias_transfer_1234_fold_4_data_type_ATAC_PE
fold_0,25M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE
fold_1,25M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_1_data_type_ATAC_PE
fold_2,25M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_2_data_type_ATAC_PE
fold_3,25M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_3_data_type_ATAC_PE
fold_4,25M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_25M/GM12878_25M_07.18.2022_bias_transfer_1234_fold_4_data_type_ATAC_PE
fold_0,5M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE
fold_1,5M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_1_data_type_ATAC_PE
fold_2,5M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_2_data_type_ATAC_PE
fold_3,5M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_3_data_type_ATAC_PE
fold_4,5M,GM12878,/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_5M/GM12878_5M_07.18.2022_bias_transfer_1234_fold_4_data_type_ATAC_PE



