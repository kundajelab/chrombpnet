#!/bin/bash
##TODO flank_size and ref_fasta files in scripts

cell_line=GM12878
date=$(date +'%m.%d.%Y')
#data_type="ATAC"
#setting=$data_type"_"$date
setting=atac_07.21.2021
cur_file_name="gm12878_testing_script.sh"

### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak
#.gz file?
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
is_filtered=True
samtools_flag=None

blacklist_region=/srv/scratch/anusri/bpnet_histone/tiledb/gm12878_dnase/onlypeaks/all_three_blacklists.bed
chrom_sizes=/srv/scratch/anusri/bpnet_histone/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/$cell_line
data_dir=$PWD/$cell_line/data
output_dir=$PWD/$cell_line/$setting


### MODEL PARAMS

gpu=5
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

## MAKE NEGATIVES BED FILE

neg_dir=$main_dir/negatives_data
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

step2_bias_json=$output_dir/bias_fit_on_signal_step2/model.0.arch
step2_bias_weights=$output_dir/bias_fit_on_signal_step2/model.0.weights


### CREATE BIGWIGS

if [[ -d $data_dir ]] ; then
    echo "skipping bigwig creation"
else
    mkdir $data_dir
    ./main_scripts/preprocess.sh $in_bam $data_dir $samtools_flag $is_filtered $data_type
    cp $PWD/$cur_file_name $data_dir
fi



### CREATE FOLDER TILEDB AND RUN TILEDB

if [[ -d $data_dir/tiledb ]] ; then
    echo "skipping tiledb"
else
    mkdir $data_dir/tiledb
    echo -e "dataset\tnegatives_peak\tidr_peak\toverlap_peak\tambig_peak\tcount_bigwig_unstranded_5p\n"$cell_line"\t"$neg_dir/bpnet.inputs.all.negatives.bed"\t"$idr_peak"\t"$overlap_peak"\t"$blacklist_region"\t"$data_dir/shifted_4_4.sorted.bam.bpnet.unstranded.bw > $data_dir/tiledb/inputs.tsv
    echo -e "overlap_peak\tbed_summit_from_last_col\nnegatives_peak\tbed_summit_from_last_col\nidr_peak\tbed_summit_from_last_col\nambig_peak\tbed_no_summit\ncount_bigwig_unstranded_5p\tbigwig" > $data_dir/tiledb/attribs.txt
    ./main_scripts/db_ingest.sh  $data_dir/tiledb/inputs.tsv $data_dir/tiledb/db $chrom_sizes $data_dir/tiledb/attribs.txt
    cp $PWD/$cur_file_name $data_dir/tiledb
fi


### STEP 1 -  TRAIN BIAS MODEL


if [[ -d $output_dir/invivo_bias_model_step1 ]] ; then
    echo "skipping step 1 - directory present "
else

    if test -z "$bias_json" 
    then
    	mkdir $output_dir/invivo_bias_model_step1

        bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "negatives_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/invivo_bias_model_step1/counts_loss_weight.txt
        counts_loss_weight_step1=`cat $output_dir/invivo_bias_model_step1/counts_loss_weight.txt`
 
   	echo -e "counts_loss_weight\t"$counts_loss_weight_step1"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/invivo_bias_model_step1/params.txt
    	params=$output_dir/invivo_bias_model_step1/params.txt
    	for fold in 0
    	do
            ./main_scripts/invivo_bias_model_step1/train.sh $fold $gpu $model_name $seed $output_dir/invivo_bias_model_step1 $params  $data_dir/tiledb/db $cell_line profile_bpnet_dnase $neg_bed_train
            ./main_scripts/invivo_bias_model_step1/predict.sh $fold $gpu $model_name $seed  $output_dir/invivo_bias_model_step1  $data_dir/tiledb/db $cell_line $chrom_sizes $neg_bed_test
            ./main_scripts/invivo_bias_model_step1/score.sh $output_dir/invivo_bias_model_step1 $model_name $fold $cell_line $seed
    	done
    	cp $PWD/$cur_file_name $output_dir/invivo_bias_model_step1
        bias_json=$output_dir/invivo_bias_model_step1/model.0.arch
        bias_weights=$output_dir/invivo_bias_model_step1/model.0.weights
    else
        echo "skipping step1 - input bias model given"
    fi

fi



### STEP 2 - FIT BIAS MODEL ON SIGNAL


if [[ -d $output_dir/bias_fit_on_signal_step2 ]] ; then
    echo "skipping step 2"
else
    mkdir $output_dir/bias_fit_on_signal_step2

    bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt
    counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt`

    counts_loss_weight_step3=$counts_loss_weight_step2

    echo -e "json_string\t"$bias_json"\nweights\t"$bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step2"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/bias_fit_on_signal_step2/params.txt
    params=$output_dir/bias_fit_on_signal_step2/params.txt
    for fold in 0
    do
        ./main_scripts/bias_fit_on_signal_step2/train.sh $fold $gpu $model_name $seed $output_dir/bias_fit_on_signal_step2 $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/bias_fit_on_signal_step2/signal_from_bias.py
        ./main_scripts/bias_fit_on_signal_step2/predict.sh $fold $gpu $model_name $seed  $output_dir/bias_fit_on_signal_step2  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/bias_fit_on_signal_step2/score.sh $output_dir/bias_fit_on_signal_step2 $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/bias_fit_on_signal_step2
fi

counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt`
counts_loss_weight_step3=$counts_loss_weight_step2


### STEP 3 - FIT BIAS AND SIGNAL MODEL

if [[ -d $output_dir/final_model_step3 ]] ; then
    echo "skipping step 3"
else
    mkdir $output_dir/final_model_step3
    echo -e "json_string\t"$step2_bias_json"\nweights\t"$step2_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3/params.txt
    params=$output_dir/final_model_step3/params.txt
    for fold in 0
    do
        ./main_scripts/final_model_step3/train.sh $fold $gpu $model_name $seed $output_dir/final_model_step3 $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/final_model_step3/profile_bpnet_dnase_with_bias.py
        ./main_scripts/final_model_step3/predict.sh $fold $gpu $model_name $seed  $output_dir/final_model_step3  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/final_model_step3/score.sh $output_dir/final_model_step3 $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3
fi


## UNPLUG MODEL

if [[ -d $output_dir/final_model_step3/unplug ]] ; then
    echo "skipping unplugging"
else
    mkdir $output_dir/final_model_step3/unplug
    unplug_bias_json=$output_dir/final_model_step3/model.0.arch
    unplug_bias_weights=$output_dir/final_model_step3/model.0.weights
    echo -e "json_string\t"$unplug_bias_json"\nweights\t"$unplug_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3/unplug/params.txt

    params=$output_dir/final_model_step3/unplug/params.txt

    for fold in 0
    do
        CUDA_VISIBLE_DEVICES=$gpu python ./main_scripts/unplug/get_model_with_bias_unplugged.py --model_params $params --outf $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 
        ./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/final_model_step3/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/unplug/score.sh $output_dir/final_model_step3/unplug $model_name $fold $cell_line $seed
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3/unplug
fi


### GET INTERPRETATIONS

if [[ -d $data_dir/$cell_line"_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $data_dir/$cell_line"_idr_split" 
    zcat $idr_peak | shuf  > $data_dir/$cell_line"_idr_split/temp.txt"
    split -l 10000 $data_dir/$cell_line"_idr_split/temp.txt" $data_dir/$cell_line"_idr_split/x"
    rm  $data_dir/$cell_line"_idr_split/temp.txt"
fi


if [[ -d $output_dir/final_model_step3/unplug/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/final_model_step3/unplug/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3/unplug/deepshap $cell_line $gpu $fold
    done

    python /srv/scratch/anusri/pipeline_script/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/final_model_step3/unplug/deepshap --target $output_dir/final_model_step3/unplug/deepshap --type 20k
    cp $PWD/$cur_file_name $output_dir/final_model_step3/unplug/deepshap

fi


### RUN MODISCO



