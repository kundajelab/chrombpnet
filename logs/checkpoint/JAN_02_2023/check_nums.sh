wc -l /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$2/$1/preds_upload/fold_0/$1_w_bias_all_regions.bed
zcat /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$2/$1/merge_folds_all_regions_may_05_24/$1_folds_merged.profile_scores_new_compressed.bed.gz | wc -l
zcat /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$2/$1/merge_folds_all_regions_may_05_24/$1_folds_merged.counts_scores_new_compressed.bed.gz | wc -l
bedtools sort -i /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$2/$1/preds_upload/fold_0/$1_w_bias_all_regions.bed | uniq | wc -l
bedtools sort -i /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$2/$1/merge_folds_all_regions_may_05_24/$1_folds_merged.counts_scores_new_compressed.bed.gz | zcat | cut -f1,2,3 | uniq | wc -l
