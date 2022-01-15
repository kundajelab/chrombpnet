kerasAC_loss_weights_bpnet  --tdb_array tiledb/db \
			    --chroms chr1 \
			    --upsample_attribute negatives_peak \
			    --label_attribute count_bigwig_unstranded_5p \
			    --num_threads 1 \
			    --task microglia \
			    --upsample_thresh 1 \
			    --flank 500
