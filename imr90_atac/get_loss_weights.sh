kerasAC_loss_weights_bpnet --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/atac \
			    --chroms chr1 \
			    --upsample_attribute overlap_peak \
			    --label_attribute count_bigwig_unstranded_5p \
			    --num_threads 1 \
			    --task IMR90 \
			    --upsample_thresh 1 \
			    --flank 500
