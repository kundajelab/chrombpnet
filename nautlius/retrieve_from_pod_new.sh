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
oak_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/$data_type_short/$dataset/$model/bias_model/bias.h5



kubectl cp $podname:/$ceph/$data_type/$dataset/$model/bias_model/bias.h5  $oak_dir

