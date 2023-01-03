OUTPREFIX=results/chrombpnet/DNASE_SE/K562/data/replicates/rep1
REFCHROMSZ=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
INTA=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/croo/K562/align/rep1/pseudorep1/K562.merged.filtered.sorted.pr1.tagAlign.gz

zcat /oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/croo/K562/align/rep1/pseudorep1/K562.merged.filtered.sorted.pr1.tagAlign.gz | sort -k1,1 -k2,2n $INTA > results/chrombpnet/DNASE_SE/K562/data/replicates/rep1.sorted.tagAlign
#gzip results/chrombpnet/DNASE_SE/K562/data/replicates/rep1.sorted.tagAlign


#bash src/helpers/preprocessing/tagaligns_to_bigwig.sh $OUTPREFIX $REFCHROMSZ results/chrombpnet/DNASE_SE/K562/data/replicates/rep1.sorted.tagAlign.gz
