#!/bin/bash
for chrom in 1 #`seq 1 22` X Y
do
    python TFscores_per_tf_calibration.py --deepSHAP ../interpret/6dil/GM12878.DNASE.bias_corrected_bpnet_tobias.fold0.deepSHAP \
	   --tf_intersection chr$chrom.GM12878.idr.peaks.overlap.tf.archetypes.bed \
	   --outf chr$chrom.GM12878.tf.deepSHAP.stats.archetypes.tsv \
	   --chrom_col 0 \
	   --summit_col 9 \
	   --peak_start_col 1 \
	   --tf_start_col 11 \
	   --tf_end_col 12 \
	   --chrom chr$chrom \
	   --vierstra_score_col 14 \
	   --tf_name_col 13
done
