#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=SURAG
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
#setting=4_$neg_shift"_shifted_"$data_type"_"$date
cur_file_name="surags_data_with_gm12878.sh"
setting=4_4_shifted_ATAC_08.12.2021

### SIGNAL INPUT

in_bigwig=/users/surag/kundajelab/scATAC-reprog/src/analysis/20210721_bias_model/data/c5_44.bw
overlap_peak=/users/surag/kundajelab/scATAC-reprog/src/analysis/20210721_bias_model/data/peaks.bed
#.gz file?
idr_peak=/users/surag/kundajelab/scATAC-reprog/src/analysis/20210721_bias_model/data/peaks.bed

blacklist_region=$PWD/data/all_three_blacklists.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/$cell_line
data_dir=$PWD/$cell_line/data
output_dir=$PWD/$cell_line/$setting


### MODEL PARAMS

gpu=0
filters=500 
n_dil_layers=8
seed=1234 
model_name=model 
neg_dir=$main_dir/negatives_data
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

## MAKE NEGATIVES BED FILE

stride=50

if [[ -d $neg_dir ]] ; then
    echo "negatives director already exists"
else
    mkdir $neg_dir
    bash main_scripts/make_gc_matched_negatives/run.sh $neg_dir $overlap_peak $ref_fasta $chrom_sizes $flank_size $stride 
fi

###  BIAS  INPUTS

bias_json=""
bias_weights=""
neg_bed_train=$neg_dir"/bpnet.inputs.train.flank300.0.negatives.bed"
neg_bed_test=$neg_dir"/bpnet.inputs.test.0.negatives.bed"

step2_bias_json=$output_dir/bias_fit_on_signal_step2_with_gm12878/model.0.arch
step2_bias_weights=$output_dir/bias_fit_on_signal_step2_with_gm12878/model.0.weights



bias_json="/srv/scratch/anusri/chrombpnet_paper/GM12878/ATAC_07.22.2021/invivo_bias_model_step1/model.0.arch"
bias_weights="/srv/scratch/anusri/chrombpnet_paper/GM12878/ATAC_07.22.2021/invivo_bias_model_step1/model.0.weights"

### STEP 2 - FIT BIAS MODEL ON SIGNAL


if [[ -d $output_dir/bias_fit_on_signal_step2_with_gm12878 ]] ; then
    echo "skipping step 2"
else
    mkdir $output_dir/bias_fit_on_signal_step2_with_gm12878

    bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/bias_fit_on_signal_step2_with_gm12878/counts_loss_weight.txt
    counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2_with_gm12878/counts_loss_weight.txt`
    counts_loss_weight_step3=$counts_loss_weight_step2

    echo -e "json_string\t"$bias_json"\nweights\t"$bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step2"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/bias_fit_on_signal_step2_with_gm12878/params.txt
    params=$output_dir/bias_fit_on_signal_step2_with_gm12878/params.txt
    for fold in 0
    do
        ./main_scripts/bias_fit_on_signal_step2/train.sh $fold $gpu $model_name $seed $output_dir/bias_fit_on_signal_step2_with_gm12878 $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/bias_fit_on_signal_step2/signal_from_bias.py
        ./main_scripts/bias_fit_on_signal_step2/predict.sh $fold $gpu $model_name $seed  $output_dir/bias_fit_on_signal_step2_with_gm12878  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/bias_fit_on_signal_step2/score.sh $output_dir/bias_fit_on_signal_step2_with_gm12878 $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/bias_fit_on_signal_step2_with_gm12878
fi


if [[ -d $data_dir/$cell_line"_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $data_dir/$cell_line"_idr_split" 
    cat $idr_peak | shuf  > $data_dir/$cell_line"_idr_split/temp.txt"
    split -l 10000 $data_dir/$cell_line"_idr_split/temp.txt" $data_dir/$cell_line"_idr_split/x"
    rm  $data_dir/$cell_line"_idr_split/temp.txt"
fi




#counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2_with_gm12878/counts_loss_weight.txt`
#counts_loss_weight_step3=$counts_loss_weight_step2


### STEP 3 - FIT BIAS AND SIGNAL MODEL

if [[ -d $output_dir/final_model_step3_with_gm12878_bias ]] ; then
    echo "skipping step 3"
else
    mkdir $output_dir/final_model_step3_with_gm12878_bias
    echo -e "json_string\t"$step2_bias_json"\nweights\t"$step2_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3_with_gm12878_bias/params.txt
    params=$output_dir/final_model_step3_with_gm12878_bias/params.txt
    for fold in 0
    do
        ./main_scripts/final_model_step3_new/train.sh $fold $gpu $model_name $seed $output_dir/final_model_step3_with_gm12878_bias $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/final_model_step3_new/profile_bpnet_dnase_with_bias.py
        ./main_scripts/final_model_step3_new/predict.sh $fold $gpu $model_name $seed  $output_dir/final_model_step3_with_gm12878_bias  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/final_model_step3_new/score.sh $output_dir/final_model_step3_with_gm12878_bias $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3_with_gm12878_bias
fi


## UNPLUG MODEL

if [[ -d $output_dir/final_model_step3_with_gm12878_bias/unplug ]] ; then
    echo "skipping unplugging"
else
    mkdir $output_dir/final_model_step3_with_gm12878_bias/unplug
    unplug_bias_json=$output_dir/final_model_step3_with_gm12878_bias/model.0.arch
    unplug_bias_weights=$output_dir/final_model_step3_with_gm12878_bias/model.0.weights
    echo -e "json_string\t"$unplug_bias_json"\nweights\t"$unplug_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3_with_gm12878_bias/unplug/params.txt

    params=$output_dir/final_model_step3_with_gm12878_bias/unplug/params.txt

    for fold in 0
    do
        CUDA_VISIBLE_DEVICES=$gpu python ./main_scripts/unplug_new/get_model_with_bias_unplugged.py --model_params $params --outf $output_dir/final_model_step3_with_gm12878_bias/unplug/$model_name.$fold.hdf5 
        ./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/final_model_step3_with_gm12878_bias/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/unplug/score.sh $output_dir/final_model_step3_with_gm12878_bias/unplug $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3_with_gm12878_bias/unplug
fi

## UNPLUG MODEL


#fold=0
#./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/final_model_step3_with_gm12878_bias/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
#./main_scripts/unplug/score.sh $output_dir/final_model_step3_with_gm12878_bias/unplug $model_name $fold $cell_line $seed

### GET INTERPRETATIONS

if [[ -d $data_dir/$cell_line"_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $data_dir/$cell_line"_idr_split" 
    cat $idr_peak | shuf  > $data_dir/$cell_line"_idr_split/temp.txt"
    split -l 10000 $data_dir/$cell_line"_idr_split/temp.txt" $data_dir/$cell_line"_idr_split/x"
    rm  $data_dir/$cell_line"_idr_split/temp.txt"
fi



if [[ -d $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_with_gm12878_bias/unplug/$model_name.$fold.hdf5 $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_with_gm12878_bias/unplug/$model_name.$fold.hdf5 $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_with_gm12878_bias/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap $cell_line $gpu $fold
    done

    python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap --target $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap --type 20k
    cp $PWD/$cur_file_name $output_dir/final_model_step3_with_gm12878_bias/unplug/deepshap

fi

modisco_sig_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
if [[ -d $modisco_sig_dir/$cell_line ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line
fi


if [[ -d $modisco_sig_dir/$cell_line/$setting"_with_gm12878_bias"/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line/$setting"_with_gm12878_bias"/
    modisco_dir_final=$modisco_sig_dir/$cell_line/$setting"_with_gm12878_bias"/
    cp  $cell_line/$setting/final_model_step3_with_gm12878_bias/unplug/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi



