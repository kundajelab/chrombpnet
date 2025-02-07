bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/modisco.inputs.bed.gz"

bedtools sort -i $bed | bedtools merge -i stdin > merged.counts.bed
bedtools sort -i $bed | bedtools merge -i stdin | bedtools intersect -a /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/k562_tf_chip_combine.bed -b stdin -wa -f 1.0 | bedtools sort -i stdin | uniq > chip_in_bed.counts.bed

#in_peaks.counts_scores_new_compressed.bed
#in_peaks.counts.interpreted_regions.bed

bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/modisco.inputs.bed.gz"

bedtools sort -i $bed | bedtools merge -i stdin > merged.profile.bed
bedtools sort -i $bed | bedtools merge -i stdin | bedtools intersect -a /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/k562_tf_chip_combine.bed -b stdin -wa -f 1.0 | bedtools sort -i stdin | uniq > chip_in_bed.profile.bed


