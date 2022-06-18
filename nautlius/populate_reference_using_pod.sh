set -beEuo pipefail

podname=populate-chrombpnet-pod
ceph=chrombpnet

echo "kubectl cp  /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits/ $podname:/$ceph/"
#kubectl cp /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits/ $podname:/$ceph/

echo "kubectl cp  /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/ $podname:/$ceph/"
#kubectl cp /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/ $podname:/$ceph/

echo "kubectl exec $podname -- mkdir -p /$ceph/chrombpnet/"
#kubectl exec $podname -- mkdir -p /$ceph/chrombpnet/

echo "kubectl cp  /mnt/lab_data2/anusri/chrombpnet/src $podname:/$ceph/chrombpnet/"
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/src  $podname:/$ceph/chrombpnet/

echo "kubectl cp  /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model.sh  $podname:/$ceph/chrombpnet/"
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model.sh  $podname:/$ceph/chrombpnet/

echo "kubectl cp  /mnt/lab_data2/anusri/chrombpnet/step4_train_bias_model.sh  $podname:/$ceph/chrombpnet/"
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step4_train_bias_model.sh  $podname:/$ceph/chrombpnet/

echo "kubectl cp  /mnt/lab_data2/anusri/chrombpnet/k562_atac_fold_0.sh  $podname:/$ceph/chrombpnet/"
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/k562_atac_fold_0.sh  $podname:/$ceph/chrombpnet/

echo "kubectl cp  /mnt/lab_data2/anusri/chrombpnet/atac_run.sh  $podname:/$ceph/chrombpnet/"
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/dnase_run.sh  $podname:/$ceph/chrombpnet/


#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run.sh  $podname:/$ceph/chrombpnet/
kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run_withbias.sh  $podname:/$ceph/chrombpnet/

#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model_dillayer.sh  $podname:/$ceph/chrombpnet/


#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model_dillayer_same_pad.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run_withbias_same_pad.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/src/training/models/chrombpnet_with_bias_model_zero_padding.py  $podname:/$ceph/chrombpnet/src/training/models/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/src/training/models/chrombpnet_with_bias_model.py  $podname:/$ceph/chrombpnet/src/training/models/


#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run_withbias_only_pos.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model_dillayer_only_pos.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/step6_train_chrombpnet_model_dillayer_same_pad_only_pos.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run_withbias_same_pad_only_pos.sh  $podname:/$ceph/chrombpnet/
#kubectl cp /mnt/lab_data2/anusri/chrombpnet/src/training/data_generators/batchgen_generator.py  $podname:/$ceph/chrombpnet/src/training/data_generators/batchgen_generator.py




#kubectl cp /mnt/lab_data2/anusri/chrombpnet/atac_run_small_negs.sh  $podname:/$ceph/chrombpnet/


#kubectl cp train.py $podname:/$ceph/chrombpnet/src/training/



