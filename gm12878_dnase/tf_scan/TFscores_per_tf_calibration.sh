#!/bin/bash
#for chrom in `seq 1 22` X Y
for chrom in 1
do
    python TFscores_per_tf_calibration.py --deepSHAP ../interpret/6dil/GM12878.DNASE.bias_corrected_bpnet_tobias.fold0.deepSHAP \
	   --tf_intersection chr$chrom.idr.peaks.overlap.tf.all.bed \
	   --outf chr$chrom.GM12878.tf.deepSHAP.stats.tsv \
	   --chrom_col 7 \
	   --summit_col 16 \
	   --peak_start_col 8 \
	   --tf_start_col 1 \
	   --tf_end_col 2 \
	   --chrom chr$chrom \
	   --vierstra_score_col 4 \
	   --tf_name_col 3
done

