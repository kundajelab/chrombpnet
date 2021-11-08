fold=0
gpu=0
model_name=model
seed=1234

bash get_predictions_in_k562_background.sh $fold $gpu $model_name $seed  $output_dir/invivo_bias_model_step1  $data_dir/tiledb/db $cell_line $chrom_sizes $neg_bed_test


