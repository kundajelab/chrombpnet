celltype=$1
oak_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$celltype/merge_folds_new/
doak_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/$celltype/merge_folds_new/
motif=$2
outf=$3
chipencid=$4
bwo=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/$celltype/data/$celltype"_unstranded.bw"
if [ $celltype = "HEPG2" ]; then
	dbwo=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$celltype/data/$celltype"_unstranded.bw"
elif [ $celltype = "K562" ]; then
	dbwo=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$celltype/data/$celltype"_unstranded.bw"
else
	dbwo=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/$celltype/data/$celltype"_unstranded.bw"
fi
echo $dbwo
python compute_imp_scores.py -ac $oak_dir/$celltype"_folds_merged.counts_scores.bw" \
	 -ap $oak_dir/$celltype"_folds_merged.profile_scores.bw" \
	 -dc $doak_dir/$celltype"_folds_merged.counts_scores.bw" \
	 -dp $doak_dir/$celltype"_folds_merged.profile_scores.bw" \
	 -cb /oak/stanford/groups/akundaje/vir/tfatlas/processed_data/$chipencid/peaks_inliers.bed.gz \
	-cc /oak/stanford/groups/akundaje/vir/tfatlas/shap/release_run_1/meanshap/$chipencid/counts_mean_shap_scores.bw \
	 -tb /mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/$celltype/$celltype"_ATAC_TOBIAS_NEW_COUNTS"/$motif/beds/$motif"_all.bed" \
	-o $outf \
	-bwo $bwo \
	-dbwo $dbwo 



