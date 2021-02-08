#for task in GM12878 H1ESC HEPG2 IMR90 K562
for task in H1ESC HEPG2 IMR90 K562
do
    genomewide_labels --task_list $task/$task.tasks.tsv \
		 --outf $task/$task.classificationlabels.SummitWithin200bpCenter.bed.gz \
		 --output_type gzip \
		 --chrom_sizes hg38.chrom.sizes \
		 --bin_stride 50 \
		 --left_flank 400 \
		 --right_flank 400 \
		 --bin_size 200 \
		 --task_threads 1 \
		 --chrom_threads 24 \
		 --allow_ambiguous \
		 --labeling_approach peak_summit_in_bin_classification \
		 --store_positives_only \
		 --split_output_by_chrom

done
