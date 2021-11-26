tdb_array=$1
chroms=$2
ambig_attribute=$3
label_attribute=$4
task=$5
flank_size=$6
out_file_name=$7

python $PWD/src/models/chrombpnet_scripts/get_loss_weights_for_bpnet.py --tdb_array $tdb_array \
                            --chroms $chroms \
                            --upsample_attribute $ambig_attribute \
                            --label_attribute $label_attribute \
                            --num_threads 1 \
                            --task  $task\
                            --upsample_thresh 1 \
                            --flank $flank_size \
                            --outf $out_file_name
                            
