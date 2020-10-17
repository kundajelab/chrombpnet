fold=0
out_dir='.'
model_prefix=tmp
CUDA_VISIBLE_DEVICES=1 kerasAC_train \
		    --seed 1234 \
		    --batch_size 50 \
		    --ref_fasta /users/annashch/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/histone \
		    --num_tasks 2 \
		    --num_inputs 3 \
		    --num_outputs 2 \
		    --architecture_spec profile_bpnet_chipseq \
		    --tdb_input_datasets seq K562,K562 K562,K562 \
		    --tdb_output_datasets K562,K562 K562,K562 \
		    --tdb_input_source_attribute seq control_count_bigwig_plus_5p,control_count_bigwig_minus_5p control_count_bigwig_plus_5p,control_count_bigwig_minus_5p \
		    --tdb_output_source_attribute count_bigwig_plus_5p,count_bigwig_minus_5p count_bigwig_plus_5p,count_bigwig_minus_5p \
		    --tdb_input_flank 3000 500,500 500,500 \
		    --tdb_output_flank 500,500 500,500 \
		    --tdb_input_min None None,None None,None \
		    --tdb_output_min None,None None,None \
		    --tdb_input_max None None,None None,None \
		    --tdb_output_max None,None None,None \
		    --tdb_input_aggregation None None,None sum,sum \
		    --tdb_input_transformation None None,None log,log \
		    --tdb_output_aggregation None,None sum,sum \
		    --tdb_output_transformation None,None log,log \
		    --tdb_partition_attribute_for_upsample overlap_peak \
		    --tdb_partition_thresh_for_upsample 1 \
		    --tdb_partition_datasets_for_upsample K562 \
		    --fold $fold \
		    --genome hg38 \
		    --num_train 100000 \
		    --num_valid 10000 \
		    --upsample_threads 24 \
		    --threads 0 \
		    --max_queue_size 100 \
		    --patience 3 \
		    --patience_lr 2 \
		    --model_prefix $out_dir/$model_prefix \
		    --model_params params.txt \
		    --upsample_ratio_list_train 1.0 \
		    --upsample_ratio_list_eval 1.0 \
		    --trackables logcount_predictions_loss loss profile_predictions_loss val_logcount_predictions_loss val_loss val_profile_predictions_loss \
		    --revcomp

