dataset=$1
model=$2
data_type=$3

scp $dataset/$model/chrombpnet_model/interpret/*.h5 anusri@sherlock.stanford.edu:/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$dataset/$model/SIGNAL/
