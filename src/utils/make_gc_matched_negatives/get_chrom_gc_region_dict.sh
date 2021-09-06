task=$1
python $PWD/main_scripts/make_gc_matched_negatives/get_chrom_gc_region_dict.py --input_bed $task/candidate.negatives.tsv --outf $task/candidate.negatives.gc.p
