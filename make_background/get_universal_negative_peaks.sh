cd /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/ENCODE_DNASE_OVERLAP_NARROWPEAK_HG38_DCC_PROCESSED_ANNASHCH
zcat *gz | bedtools sort -i - | bedtools merge -i - | bedtools subtract -a ~/hg38.chrom.sizes.bed -b - > /srv/scratch/annashch/chrombpnet/make_background/universal_negatives.bed

