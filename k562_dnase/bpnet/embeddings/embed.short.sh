for fold in 0
do
    CUDA_VISIBLE_DEVICES=0 compute_nn_embeddings \
			--input_bed_file enhancers.bed \
			--model_hdf5 /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.$fold.hdf5 \
			--ref_fasta /mnt/data/male.hg19.fa \
			--flank 132 \
			--center_on_bed_interval \
			--output_npz_file jesse_k562_fold0_dil5.npz \
			--embedding_layer_name 5_dilconv \
			--input_layer_name sequence \
			--threads 15 \
			--batch_size 25 \
			--global_pool_on_position
done
