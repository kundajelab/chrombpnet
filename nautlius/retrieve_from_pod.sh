#set -beEuo pipefail


dataset=$1
model=$2
data_type=$3

podname=populate-chrombpnet-pod
ceph=chrombpnet


mkdir /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/$dataset/nautilus_runs_jun16/
echo "kubectl cp $podname:/$ceph/$data_type/$dataset/$model  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/$dataset/nautilus_runs_jun16/"
kubectl cp $podname:/$ceph/$data_type/$dataset/$model  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/$dataset/nautilus_runs_jun16/




#mkdir /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$dataset/nautilus_runs_may18/
#echo "kubectl cp $podname:/$ceph/chrombpnet/DNASE_PE/$dataset/  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$dataset/nautilus_runs_apr12/"
#echo "kubectl cp $podname:/$ceph/chrombpnet/DNASE_PE/$dataset/$dataset_04.09.2022_bias_128_4_1234_0.4_fold_0  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$dataset/nautilus_runs_apr12/"
#kubectl cp $podname:/$ceph/DNASE_PE/$dataset/$dataset"_05.13.2022_bias_128_4_1234_0.5_fold_0"  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$dataset/nautilus_runs_may18/$dataset"_05.13.2022_bias_128_4_1234_0.5_fold_0"


#kubectl cp $podname:/$ceph/ATAC_PE/$dataset/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_2/  /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/$dataset/nautilus_runs_apr12/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_2/

#model=K562_05.09.2022_bias_128_4_1234_1.0_fold_0
#model=K562_05.10.2022_bias_128_4_1234_0.8_fold_0
#model=H1ESC_05.09.2022_bias_128_4_1234_0.8_fold_0
#model=HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
#model=H1ESC_04.09.2022_bias_128_4_1234_0.5_fold_0
#model=GM12878_05.10.2022_bias_128_4_1234_0.8_fold_0
#model=HEPG2_06.06.2022_bias_128_8_1234_0.8_fold_0



#mkdir /oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/DNASE_PE/$dataset/$model/
#mkdir /oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/DNASE_PE/$dataset/$model/SIGNAL/
#oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/DNASE_PE/$dataset/$model/SIGNAL

#kubectl cp $podname:/$ceph/DNASE_PE/$dataset/$model/bias_model/interpret/$dataset.counts_scores.h5  $oak_dir/$dataset.counts_scores.h5
#wait
#kubectl cp $podname:/$ceph/DNASE_PE/$dataset/$model/bias_model/interpret/$dataset.profile_scores.h5  $oak_dir/$dataset.profile_scores.h5

#kubectl cp $podname:/$ceph/DNASE_PE/$dataset/$model/chrombpnet_model/interpret/$dataset.counts_scores.h5  $oak_dir/$dataset.counts_scores.h5
#wait

#kubectl cp $podname:/$ceph/DNASE_PE/$dataset/$model/chrombpnet_model/interpret/$dataset.profile_scores.h5  $oak_dir/$dataset.profile_scores.h5


