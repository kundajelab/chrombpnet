#bw1=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/interpret/merged.K562.profile.bw
#bw1=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/background_interpret/K562.profile.bw
#bed1=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/background_interpret/K562.interpreted_regions_v2.bed

#head $bed1
#K562_wo_bias.bw
#full_K562.profile.bigwig 

#bw2=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_1_data_type_ATAC_PE/chrombpnet_model/interpret/full_K562.profile.bigwig
#bw2=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_1_data_type_ATAC_PE/background_interpret/K562.profile.bw
#bed2=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_1_data_type_ATAC_PE/background_interpret/K562.interpreted_regions_v2.bed

#head $bed2

#bw3=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_2_data_type_ATAC_PE/chrombpnet_model/interpret/full_K562.profile.bigwig
#bw3=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_2_data_type_ATAC_PE/background_interpret/K562.profile.bw
#bed3=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_2_data_type_ATAC_PE/background_interpret/K562.interpreted_regions_v2.bed

#head $bed3

#bw4=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_3_data_type_ATAC_PE/chrombpnet_model/interpret/full_K562.profile.bigwig
#bw4=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_3_data_type_ATAC_PE/background_interpret/K562.profile.bw
#bed4=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_3_data_type_ATAC_PE/background_interpret/K562.interpreted_regions_v2.bed

#head $bed4

#bw5=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_4_data_type_ATAC_PE/chrombpnet_model/interpret/full_K562.profile.bigwig

#bw5=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_4_data_type_ATAC_PE/background_interpret/K562.profile.bw
#bed5=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/K562_07.07.2022_bias_128_4_2356_0.5_fold_4_data_type_ATAC_PE/background_interpret/K562.interpreted_regions_v2.bed

#head $bed5

#output=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/merge_folds/K562.profile.background.bw

bw1=$1
bw2=$2
bw3=$3
bw4=$4
bw5=$5
output=$6

wiggletools mean $bw1 $bw2 $bw3 $bw4 $bw5 | wigToBigWig stdin /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes $output



