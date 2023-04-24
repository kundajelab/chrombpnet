setting=$1
fold=$2
mkdir  /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/HEPG2_DNASE_PE/$fold/
kubectl cp populate-chrombpnet-pod:/chrombpnet/DNASE_PE/HEPG2/$setting/bias_model/bias.h5 /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/$setting/bias_model/bias.h5
kubectl cp populate-chrombpnet-pod:/chrombpnet/DNASE_PE/HEPG2/$setting/bias_model/bias.h5 /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/HEPG2_DNASE_PE/$fold/bias.h5


