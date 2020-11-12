for fold in 0
do
    CUDA_VISIBLE_DEVICES=0 compute_nn_embeddings \
			--input_bed_file enhancers.bed \
			--model_hdf5 /srv/scratch/anusri/anna_k562_dnase/k562.dnase.with.bias.unplugged.$fold.hdf5 \
			--ref_fasta /mnt/data/male.hg19.fa \
			--flank 132 \
			--output_npz_file k562_enhancer_embeddings_$fold.5_dilconv.npz \
			--embedding_layer_name 5_dilconv \
			--input_layer_name sequence \
                        --center_on_bed_interval True \
			--threads 20 \
			--batch_size 50 \
			--global_pool_on_position
done
