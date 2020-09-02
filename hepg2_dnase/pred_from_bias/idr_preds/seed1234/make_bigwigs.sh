kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds/tobias.bias.preds.hepg2.0.predictions \
			--out_prefix hepg2_dnase_bpnet \
			--chrom_sizes_file hg38.chrom.sizes &
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds/6mer.bias.preds.hepg2.0.predictions \
			--out_prefix hepg2_dnase_tobias \
			--chrom_sizes_file hg38.chrom.sizes &
kerasAC_bigwigs_from_io --predictions_hdf /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds/bpnet.bias.preds.hepg2.0.predictions \
			--out_prefix hepg2_dnase_6mer \
			--chrom_sizes_file hg38.chrom.sizes &


