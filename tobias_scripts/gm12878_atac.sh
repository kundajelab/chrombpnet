cell_line=GM12878
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
genome_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
overlap=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak

blacklist=$PWD/data/all_three_blacklists.bed
outdir=$PWD/tobias_scripts/$cell_line

bedtools slop -i $overlap -g $PWD/data/hg38.chrom.sizes -b 1057 > $PWD/tobias_scripts/$cell_line/overlap_flank1057.bed
sort -k1,1 -k2,2n $PWD/tobias_scripts/$cell_line/overlap_flank1057.bed > $PWD/tobias_scripts/$cell_line/sorted.bed
bedtools merge -i $PWD/tobias_scripts/$cell_line/sorted.bed > $PWD/tobias_scripts/$cell_line/merged.bed
merged_overlap=$PWD/tobias_scripts/$cell_line/merged.bed

bash $PWD/tobias_scripts/tobias.sh $in_bam $genome_fasta $merged_overlap $blacklist $outdir
