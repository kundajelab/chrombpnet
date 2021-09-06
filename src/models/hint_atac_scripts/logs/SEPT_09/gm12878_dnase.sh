cell_line=GM12878/DNASE
in_bam=/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam
overlap=/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz

bedtools slop -i $overlap -g $PWD/data/hg38.chrom.sizes -b 1057 > $PWD/hint_atac_scripts/$cell_line/overlap_flank1057_dnase.bed
sort -k1,1 -k2,2n $PWD/hint_atac_scripts/$cell_line/overlap_flank1057_dnase.bed > $PWD/hint_atac_scripts/$cell_line/sorted_dnase.bed
bedtools merge -i $PWD/hint_atac_scripts/$cell_line/sorted_dnase.bed > $PWD/hint_atac_scripts/$cell_line/merged_dnase.bed
merged_overlap=$PWD/hint_atac_scripts/$cell_line/merged_dnase.bed

bash $PWD/hint_atac_scripts/hint_dnase.sh $cell_line $in_bam $merged_overlap
