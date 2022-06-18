encode_id="ENCSR108PUO"
chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR108PUO/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR108PUO/chrombpnet_model/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR108PUO/chrombpnet_model/bias_model_scaled.h5
gpu=0
regions=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/summit_ccres.tsv
output_dir=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/pgp1/
celline=$encode_id

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


