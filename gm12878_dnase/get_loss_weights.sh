#counts_loss_weight 64.0
kerasAC_loss_weights_bpnet  --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/dnase \
			    --chroms chr1 \
			    --upsample_attribute idr_peak \
			    --label_attribute count_bigwig_unstranded_5p \
			    --num_threads 1 \
			    --task ENCSR000EMT \
			    --upsample_thresh 1 \
			    --flank 500


