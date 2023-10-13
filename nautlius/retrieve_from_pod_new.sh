#set -beEuo pipefail


dataset=$1
model=$2
data_type=$3
data_type_short=$4

podname=populate-chrombpnet-pod
ceph=chrombpnet

#mkdir /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/
#mkdir /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/$model
#mkdir /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/$model/bias_model
#oak_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/$model/bias_model/bias.h5
oak_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/$model/bias_model/
mkdir $oak_dir



echo "kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/  $oak_dir"



filen=bias_metrics.json
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias.args.json
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias_data_params.tsv
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias.h5
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias.log
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias.log.batch
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias_model_params.tsv
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=bias.params.json
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=filtered.bias_nonpeaks.bed
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen

filen=train_bias_model.log
kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/$filen  $oak_dir/$filen







