model=/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_nov_08/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/unplug/model.0.hdf5
python snp_scoring.py -i /mnt/lab_data2/anusri/variant_effect_prediction_example/subsample_test.csv -g /mnt/data/male.hg19.fa -m $model -o $PWD -bs 64 --debug_mode_on 0
#python snp_scoring.py -i subsample_test_new_3.csv -g /mnt/data/male.hg19.fa -m $model -o $PWD -bs 64 --debug_mode_on 1

