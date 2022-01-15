#generate candidate negative regions
cat overlap.optimal_peak.narrowPeak /mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed | cut -f1,2,3 > non.negatives.bed
bedtools subtract -a ~/hg38.chrom.sizes.bed -b non.negatives.bed > candidate.negatives.bed


prefix=/users/annashch/kerasAC/kerasAC/helpers/generate_training_intervals

python $prefix/bed_to_model_input_windows.py --bed ../data/microglia.overlap.optimal.narrowPeak --chrom_sizes ~/hg38.chrom.sizes --outf positive.expanded.bed
python $prefix/bed_to_model_input_windows.py --bed candidate.negatives.bed --chrom_sizes ~/hg38.chrom.sizes --use_center --dont_expand_to_full_peak --outf candidate.negatives.expanded.bed

python $prefix/get_gc_content.py --input_bed positive.expanded.bed --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --out_pickle gc.positives.pkl
python $prefix/get_gc_content.py --input_bed candidate.negatives.expanded.bed --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --out_pickle gc.negatives.pkl
python $prefix/get_gc_matched_negatives.py --gc_content_pickle_positives gc.positives.pkl --gc_content_pickle_negatives gc.negatives.pkl --outf microglia.gc.matched.negatives.bed



