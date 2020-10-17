for fold in `seq 0 4`
do
    CUDA_VISIBLE_DEVICES=1 compute_nn_embeddings \
			--input_bed_file /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/embeddings/idr.optimal_peak.narrowPeak.gz \
			--model_hdf5 /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.$fold.hdf5 \
			--ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
			--center_on_summit \
			--flank 673 \
			--output_npz_file k562_dnase_unplugged_embeddings_$fold.5_dilconv.npz \
			--embedding_layer_name 5_dilconv \
			--input_layer_name sequence \
			--threads 15 \
			--batch_size 25 \
			--global_pool_on_position
done
