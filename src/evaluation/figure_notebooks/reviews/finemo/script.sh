shuf /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/filtered.peaks.bed | head -n 50000 > random_1000_peaks.bed

file1=/mnt/lab_data2/anusri/finemo_gpu/output/572M/hits_unique.tsv
file3=random_1000_peaks.bed


tail -n +2 "$file1" |
awk 'BEGIN {OFS="\t"} {print $1, $4, $5, $2, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' |
bedtools intersect -a stdin -b "$file3" -wa -f 1.0 > intersected.bed


bedtools sort -i random_1000_peaks.bed | bedtools merge -i stdin  > merged_peaks.bed
bedtools sort -i intersected.bed | bedtools merge -i stdin > merged_intersected.bed

bedtools slop -i merged_intersected.bed -g /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes -b 5 > slop_intersect.bed


cat intersected.bed | 
awk 'BEGIN {OFS="\t"} {print $1, $4, $5, $2, $3, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' |
bedtools sort -i stdin |
bedtools merge -i stdin > hits_inpeaks.bed

bedtools subtract -a merged_peaks.bed -b slop_intersect.bed > outside_hits.bed


rm intersected.bed
rm slop_intersect.bed
rm merged_intersected.bed
rm merged_peaks.bed

