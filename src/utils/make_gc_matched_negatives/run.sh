#!/bin/bash

foreground_bed=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
exclude_bed=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
flank_size=1057
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
output_dir=output/
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
genomewide_gc=/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38.tsv

echo "get gc content of the positive sequences" 
python $PWD/src/utils/make_gc_matched_negatives/get_gc_content.py \
        --input_bed $foreground_bed \
        --ref_fasta $ref_fasta \
        --out_prefix $output_dir/foreground.gc.bed \
        --flank_size $flank_size \

echo "get candidate negative bed" 
bedtools intersect -v -a $genomewide_gc -b $exclude_bed  > $output_dir/candidate.negatives.tsv

echo "find regions in candidate negative bed that gc-match wth foreground" 
bash $PWD/src/utils/make_gc_matched_negatives/get_gc_matched_negatives.py \
        --candidate_negatives $output_dir/candidate.negatives.tsv \
        --foreground_gc_bed $output_dir/foreground.gc.bed \
        --out_prefix $output_dir/negatives.bed

