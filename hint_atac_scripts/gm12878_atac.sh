cell_line=GM12878
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
overlap=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak

bedtools slop -i $overlap -g $PWD/data/hg38.chrom.sizes -b 1057 > $PWD/hint_atac_scripts/$cell_line/overlap_flank1057.bed
sort -k1,1 -k2,2n $PWD/hint_atac_scripts/$cell_line/overlap_flank1057.bed > $PWD/hint_atac_scripts/$cell_line/sorted.bed
bedtools merge -i $PWD/hint_atac_scripts/$cell_line/sorted.bed > $PWD/hint_atac_scripts/$cell_line/merged.bed
merged_overlap=$PWD/hint_atac_scripts/$cell_line/merged.bed

bash hint_atac.sh $cell_line $in_bam $merged_overlap
