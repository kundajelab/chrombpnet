
cellty=GM12878
#cp jan_3_2024/mean_modisco_plots/ATAC/$cellty/chip_in_bed.bed oct_24_2024/ATAC/$cellty/


counts=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$cellty/merge_folds_new_may_05_24/counts/modisco_counts.h5
profile=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$cellty/merge_folds_new_may_05_24/profile/modisco_profile.h5

ocounts=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$cellty/merge_folds_new_may_05_24/counts/modisco_old_format.h5
oprofile=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/$cellty/merge_folds_new_may_05_24/profile/modisco_old_format.h5

modisco convert-backward -i $counts -o $ocounts
modisco convert-backward -i $profile -o $oprofile



