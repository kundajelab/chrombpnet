set -beEuo pipefail

experiment=$1
data_type=DNASE_SE
#data_type=ATAC_PE

#oak_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet_feb_04/$data_type/
oak_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/


podname=populate-chrombpnet-pod
ceph=chrombpnet

echo "kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment"
#kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment


echo "kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment/data"
#kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment/data


#echo "kubectl cp $oak_dir/$experiment/data/$experiment""_unstranded.bw $podname:/$ceph/$data_type/$experiment/data/"
#kubectl cp $oak_dir/$experiment/data/$experiment"_unstranded.bw" $podname:/$ceph/$data_type/$experiment/data/


#echo "kubectl cp $oak_dir/$experiment/data/peaks_no_blacklist.bed $podname:/$ceph/$data_type/$experiment/data/"
#kubectl cp $oak_dir/$experiment/data/peaks_no_blacklist.bed $podname:/$ceph/$data_type/$experiment/data/


#echo "kubectl cp $oak_dir/$experiment/data/30K.subsample.overlap.bed $podname:/$ceph/$data_type/$experiment/data/"
#kubectl cp $oak_dir/$experiment/data/30K.subsample.overlap.bed $podname:/$ceph/$data_type/$experiment/data/



#echo "kubectl cp $oak_dir/$experiment/negatives_data/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/"
#kubectl cp $oak_dir/$experiment/negatives_data/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_30_stride.bed
#kubectl cp $oak_dir/$experiment/negatives_data/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit.bed

#kubectl cp $oak_dir/$experiment/testing_1to1_negatives_data $podname:/$ceph/$data_type/$experiment/testing_1to1_negatives_data/
kubectl cp $oak_dir/$experiment/negatives_data_1/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_1.bed
kubectl cp $oak_dir/$experiment/negatives_data_2/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_2.bed
kubectl cp $oak_dir/$experiment/negatives_data_3/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_3.bed
kubectl cp $oak_dir/$experiment/negatives_data_4/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_4.bed










