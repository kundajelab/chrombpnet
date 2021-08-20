#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=GM12878
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
#setting=4_$neg_shift"_shifted_"$data_type"_"$date
cur_file_name="withhepg2bias_gm12878_atac_script.sh"
setting=ATAC_07.22.2021
transfer_cell_line=HEPG2

### SIGNAL INPUT

blacklist_region=$PWD/data/all_three_blacklists.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/$cell_line
data_dir=$PWD/$cell_line/data
output_dir=$PWD/$cell_line/$setting


### MODEL PARAMS

gpu=1
filters=500 
n_dil_layers=8
seed=1234 
model_name=model 
neg_dir=$main_dir/neg_data
flank_size=1057

## CREATE DIRS

if [[ -d $main_dir ]] ; then
    echo "main director already exists"
else
    mkdir $main_dir
fi

if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir $output_dir
fi



### STEP 2 - FIT BIAS MODEL ON SIGNAL
bias_json=$PWD/$transfer_cell_line/$setting/invivo_bias_model_step1/model.0.arch
bias_weights=$PWD/$transfer_cell_line/$setting/invivo_bias_model_step1/model.0.weights

if [[ -d $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2" ]] ; then
    echo "skipping step 2"
else
    mkdir $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"

    bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/counts_loss_weight.txt
    counts_loss_weight_step2=`cat $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/counts_loss_weight.txt`
    counts_loss_weight_step3=$counts_loss_weight_step2

    echo -e "json_string\t"$bias_json"\nweights\t"$bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step2"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/params.txt
    params=$output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/params.txt
    for fold in 0
    do
        ./main_scripts/bias_fit_on_signal_step2/train.sh $fold $gpu $model_name $seed $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2" $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/bias_fit_on_signal_step2/signal_from_bias.py
        ./main_scripts/bias_fit_on_signal_step2/predict.sh $fold $gpu $model_name $seed  $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/bias_fit_on_signal_step2/score.sh $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2" $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"
fi

counts_loss_weight_step2=`cat $output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/counts_loss_weight.txt`
counts_loss_weight_step3=$counts_loss_weight_step2

### STEP 3 - FIT BIAS AND SIGNAL MODEL

step2_bias_json=$output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/model.0.arch
step2_bias_weights=$output_dir/$transfer_cell_line"_bias_fit_on_signal_step2"/model.0.weights

if [[ -d $output_dir/with_$transfer_cell_line"_bias_final_model" ]] ; then
    echo "skipping step 3"
else
    mkdir $output_dir/with_$transfer_cell_line"_bias_final_model"
    echo -e "json_string\t"$step2_bias_json"\nweights\t"$step2_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/with_$transfer_cell_line"_bias_final_model"/params.txt
    params=$output_dir/with_$transfer_cell_line"_bias_final_model"/params.txt
    for fold in 0
    do
        ./main_scripts/final_model_step3_new/train.sh $fold $gpu $model_name $seed $output_dir/with_$transfer_cell_line"_bias_final_model" $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/final_model_step3_new/profile_bpnet_dnase_with_bias.py
        ./main_scripts/final_model_step3_new/predict.sh $fold $gpu $model_name $seed  $output_dir/with_$transfer_cell_line"_bias_final_model"  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/final_model_step3_new/score.sh $output_dir/with_$transfer_cell_line"_bias_final_model" $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/with_$transfer_cell_line"_bias_final_model"
fi


## UNPLUG MODEL

if [[ -d $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug ]] ; then
    echo "skipping unplugging"
else
    mkdir $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug
    unplug_bias_json=$output_dir/with_$transfer_cell_line"_bias_final_model"/model.0.arch
    unplug_bias_weights=$output_dir/with_$transfer_cell_line"_bias_final_model"/model.0.weights
    echo -e "json_string\t"$unplug_bias_json"\nweights\t"$unplug_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/params.txt

    params=$output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/params.txt

    for fold in 0
    do
        CUDA_VISIBLE_DEVICES=$gpu python ./main_scripts/unplug_new/get_model_with_bias_unplugged.py --model_params $params --outf $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/$model_name.$fold.hdf5 
        ./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/unplug/score.sh $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug
fi


#fold=0
#./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
#./main_scripts/unplug/score.sh $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug $model_name $fold $cell_line $seed

### GET INTERPRETATIONS

if [[ -d $data_dir/$cell_line"_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $data_dir/$cell_line"_idr_split" 
    zcat $idr_peak | shuf  > $data_dir/$cell_line"_idr_split/temp.txt"
    split -l 10000 $data_dir/$cell_line"_idr_split/temp.txt" $data_dir/$cell_line"_idr_split/x"
    rm  $data_dir/$cell_line"_idr_split/temp.txt"
fi


if [[ -d $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret.sh $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/$model_name.$fold.hdf5 $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/$model_name.$fold.hdf5 $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap $cell_line $gpu $fold
    done

    python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap --target $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap --type 20k
    cp $PWD/$cur_file_name $output_dir/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap

fi

modisco_sig_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
if [[ -d $modisco_sig_dir/$cell_line ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line
fi


if [[ -d $modisco_sig_dir/$cell_line/$setting"_with_"$transfer_cell_line"_bias"/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line/$setting"_with_"$transfer_cell_line"_bias"/
    modisco_dir_final=$modisco_sig_dir/$cell_line/$setting"_with_"$transfer_cell_line"_bias"/
    cp  $cell_line/$setting/with_$transfer_cell_line"_bias_final_model"/unplug/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi



