encode_id="ENCSR330JXM"
chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/$encode_id/chrombpnet_model_encsr146kfx_bias/chrombpnet_wo_bias.h5
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/$encode_id/chrombpnet_model_encsr146kfx_bias/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/$encode_id/chrombpnet_model_encsr146kfx_bias/bias_model_scaled.h5
gpu=3
regions=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/summit_ccres.tsv
output_dir=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/brain_90d/
celline=$encode_id

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


