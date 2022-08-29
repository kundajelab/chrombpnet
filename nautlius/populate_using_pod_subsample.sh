set -beEuo pipefail



depth=$1

experiment=GM12878_$depth
oak_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data_$depth/
data_type=ATAC_PE

podname=populate-chrombpnet-pod
ceph=chrombpnet


echo "kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment"
kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment


echo "kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment/data"
kubectl exec $podname -- mkdir -p /$ceph/$data_type/$experiment/data


echo "kubectl cp $oak_dir/"GM12878_unstranded.bw" $podname:/$ceph/$data_type/$experiment/data/"
kubectl cp $oak_dir/"GM12878_unstranded.bw" $podname:/$ceph/$data_type/$experiment/data/


echo "kubectl cp $oak_dir/peaks_no_blacklist.bed $podname:/$ceph/$data_type/$experiment/data/"
kubectl cp $oak_dir/peaks_no_blacklist.bed $podname:/$ceph/$data_type/$experiment/data/


echo "kubectl cp /mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed $podname:/$ceph/$data_type/$experiment/data/"
kubectl cp /mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed $podname:/$ceph/$data_type/$experiment/data/


kubectl cp $oak_dir/negatives_data/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit.bed
kubectl cp $oak_dir/negatives_data_1/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_1.bed
kubectl cp $oak_dir/negatives_data_2/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_2.bed
kubectl cp $oak_dir/negatives_data_3/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_3.bed
kubectl cp $oak_dir/negatives_data_4/negatives_with_summit.bed $podname:/$ceph/$data_type/$experiment/data/negatives_with_summit_4.bed










