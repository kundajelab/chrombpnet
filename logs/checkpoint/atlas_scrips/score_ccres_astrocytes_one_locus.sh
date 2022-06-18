encode_id="ENCSR146KFX"
chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR146KFX/chrombppnet_model_encsr283tme_bias/chrombpnet_wo_bias.h5
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR146KFX/chrombppnet_model_encsr283tme_bias/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR146KFX/chrombppnet_model_encsr283tme_bias/bias_model_scaled.h5
#chrombpnet_nb=/oak/stanford/groups/akundaje/projects/neuro-variants/outs/1_31_2022_adpd_model_training/full_models/cluster13_fold0/chrombpnet_wo_bias.h5
#chrombpnet=/oak/stanford/groups/akundaje/projects/neuro-variants/outs/1_31_2022_adpd_model_training/full_models/cluster13_fold0/chrombpnet.h5
#bias=/oak/stanford/groups/akundaje/projects/neuro-variants/outs/1_31_2022_adpd_model_training/full_models/cluster13_fold0/bias_model_scaled.h5
gpu=7
regions=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/astrocytes_one_locus_peaks.tsv
#output_dir=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/astrocytes_one_locus/single_cell/
output_dir=/mnt/lab_data2/anusri/chromatin_atlas_dnase/brain_organoid_dev/astrocytes_one_locus/
celline=$encode_id

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


