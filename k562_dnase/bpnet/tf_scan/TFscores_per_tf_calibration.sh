#!/bin/bash
#for chrom in `seq 1 22` X Y
for chrom in 1
do
    python TFscores_per_tf_calibration.py --deepSHAP /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/interpret/6dil/K562.DNASE.bias_corrected_bpnet_tobias.unplugged.fold0.deepSHAP \
	   --tf_intersection chr$chrom.K562.idr.peaks.overlap.tf.all.bed \
	   --outf chr$chrom.K562.tf.deepSHAP.stats.tsv \
	   --chrom_col 0 \
	   --summit_col 9 \
	   --peak_start_col 1 \
	   --tf_start_col 11 \
	   --tf_end_col 12 \
	   --chrom chr$chrom \
	   --vierstra_score_col 14 \
	   --tf_name_col 13
done
