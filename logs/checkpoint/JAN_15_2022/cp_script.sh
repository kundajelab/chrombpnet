data_type=DNASE_SE
cell_line=$1
#setting=ATAC_PE_12.30.2021
#setting=ATAC_PE_withk562bias_01.03.2022
setting=DNASE_SE_withk562bias_01.04.2022
oak_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/$data_type/$cell_line/
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
output_dir=$main_dir/$setting
if [[ -d $oak_dir/$setting/SIGNAL ]]; then
    echo "dir exists"
else
    mkdir $oak_dir/$setting/
    mkdir $oak_dir/$setting/SIGNAL
fi

if [[ -f $oak_dir/$setting/SIGNAL/$cell_line.counts_scores.h5  && -f $oak_dir/$setting/SIGNAL/$cell_line.profile_scores.h5 ]] ; then
    echo "files present"
else
    echo "copying counts and profile scores to oak"
    cp $output_dir/chrombpnet_model/interpret/$cell_line.counts_scores.h5 $oak_dir/$setting/SIGNAL/
    cp $output_dir/chrombpnet_model/interpret/$cell_line.profile_scores.h5 $oak_dir/$setting/SIGNAL/
fi
