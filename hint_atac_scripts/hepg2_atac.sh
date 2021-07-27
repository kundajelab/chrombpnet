cell_line=HEPG2
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/HEPG2.atac.filt.merged.bam
overlap=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak

bedtools slop -i $overlap -g $PWD/data/hg38.chrom.sizes -b 1057 > $PWD/hint_atac_scripts/$cell_line/overlap_flank1057.bed
sort -k1,1 -k2,2n $PWD/hint_atac_scripts/$cell_line/overlap_flank1057.bed > $PWD/hint_atac_scripts/$cell_line/sorted.bed
bedtools merge -i $PWD/hint_atac_scripts/$cell_line/sorted.bed > $PWD/hint_atac_scripts/$cell_line/merged.bed
merged_overlap=$PWD/hint_atac_scripts/$cell_line/merged.bed

bash hint_atac.sh $cell_line $in_bam $merged_overlap
