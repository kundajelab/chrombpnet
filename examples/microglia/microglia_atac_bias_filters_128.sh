#!/bin/bash
#g#TODO flank_size and ref_fasta files in scripts 


cell_line=microglia
data_type="ATAC"
neg_shift=4
filters=128
n_dil_layers=4

date=$(date +'%m.%d.%Y')
setting=4_$neg_shift"_shifted_"$data_type"_"$date"_bias_filters_"$filters
cur_file_name="microglia_atac_bias_filters_"$filters".sh"
setting=4_4_shifted_ATAC_10.11.2021_bias_filters_128

### SIGNAL INPUT

#in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
#overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak
#idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
#is_filtered=True
#samtools_flag=None

#blacklist_region=/data/GRch38_unified_blacklist.bed
chrom_sizes=~/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
mkdir -p $main_dir
data_dir=$PWD
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting
mkdir -p $output_dir

### MODEL PARAMS

gpu=0
seed=1234 
model_name=microglia 
neg_dir=$PWD/negatives
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

#stride=50
#if [[ -d $neg_dir ]] ; then
#    echo "negatives director already exists"
#else
    mkdir $neg_dir
#    bash $PWD/src/utils/make_gc_matched_negatives/run.sh $neg_dir $overlap_peak $ref_fasta $chrom_sizes $flank_size $stride $blacklist_region 
#fi


###  BIAS  INPUTS

bias_json=""
bias_weights=""
neg_bed_train=$neg_dir"/train.fold0.negatives.bed"
neg_bed_test=$neg_dir"/test.fold0.negatives.bed"


### CREATE BIGWIGS

#if [[ -d $data_dir ]] ; then
#    echo "skipping bigwig creation"
#else
#    mkdir $data_dir
#    bash $PWD/src/utils/preprocess.sh $in_bam $data_dir $samtools_flag $is_filtered $data_type $neg_shift
#    cp $PWD/$cur_file_name $data_dir
#fi


### CREATE FOLDER TILEDB AND RUN TILEDB

#if [[ -d $data_dir/tiledb ]] ; then
#    echo "skipping tiledb"
#else
#    mkdir $data_dir/tiledb
#    echo -e "dataset\tnegatives_peak\tidr_peak\toverlap_peak\tambig_peak\tcount_bigwig_unstranded_5p\n"$cell_line"\t"$neg_dir/bpnet.inputs.all.negatives.bed"\t"$idr_peak"\t"$overlap_peak#"\t"$blacklist_region"\t"$data_dir/shifted.sorted.bam.bpnet.unstranded.bw > $data_dir/tiledb/inputs.tsv
#    echo -e "overlap_peak\tbed_summit_from_last_col\nnegatives_peak\tbed_summit_from_last_col\nidr_peak\tbed_summit_from_last_col\nambig_peak\tbed_no_summit\ncount_bigwig_unstranded_5p\tbigwig" > $data_dir/tiledb/attribs.txt
#    bash $PWD/src/utils/db_ingest.sh  $data_dir/tiledb/inputs.tsv $data_dir/tiledb/db $chrom_sizes $data_dir/tiledb/attribs.txt
#    cp $PWD/$cur_file_name $data_dir/tiledb
#fi



### STEP 1 -  TRAIN BIAS MODEL


if [[ -d $output_dir/invivo_bias_model_step1 ]] ; then
    echo "skipping step 1 - directory present "
else

    if test -z "$bias_json" 
    then
        mkdir $output_dir/invivo_bias_model_step1

        #bash $PWD/src/models/chrombpnet_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "negatives_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/invivo_bias_model_step1/counts_loss_weight.txt
        #counts_loss_weight_step1=`cat $output_dir/invivo_bias_model_step1/counts_loss_weight.txt`
 
        #echo -e "counts_loss_weight\t"$counts_loss_weight_step1"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/invivo_bias_model_step1/params.txt
	echo -e "counts_loss_weight\t2.1\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/invivo_bias_model_step1/params.txt
        params=$output_dir/invivo_bias_model_step1/params.txt
        for fold in 0
        do
            min_logcount=2.3
            max_logcount=11.5
            bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/train.sh $fold $gpu $model_name $seed $output_dir/invivo_bias_model_step1 $params  $data_dir/tiledb/db $cell_line $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/bpnet_model.py $neg_bed_train $min_logcount $max_logcount
            bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/predict.sh $fold $gpu $model_name $seed  $output_dir/invivo_bias_model_step1  $data_dir/tiledb/db $cell_line $chrom_sizes $neg_bed_test
            bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/score.sh $output_dir/invivo_bias_model_step1 $model_name $fold $cell_line $seed $min_logcount $max_logcount
        done
        cp $PWD/$cur_file_name $output_dir/invivo_bias_model_step1
        bias_json=$output_dir/invivo_bias_model_step1/microglia.0.arch
        bias_weights=$output_dir/invivo_bias_model_step1/microglia.0.weights
    else
        echo "skipping step1 - input bias model given"
    fi

fi

fold=0
#bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/predict_inpeak.sh $fold $gpu $model_name $seed  $output_dir/invivo_bias_model_step1  $data_dir/tiledb/db $cell_line $chrom_sizes
#bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/score_inpeak.sh $output_dir/invivo_bias_model_step1 $model_name $fold $cell_line $seed
#bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/score.sh $output_dir/invivo_bias_model_step1 $model_name $fold $cell_line $seed


filters=500 
n_dil_layers=8

step2_bias_json=$output_dir/invivo_bias_model_step1/microglia.0.arch
step2_bias_weights=$output_dir/invivo_bias_model_step1/microglia.0.weights
### STEP 3 - FIT BIAS AND SIGNAL MODEL

if [[ -d $output_dir/final_model_step3 ]] ; then
    echo "skipping step 3"
else
   mkdir $output_dir/final_model_step3
   bash $PWD/src/models/chrombpnet_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/final_model_step3/counts_loss_weight.txt
    counts_loss_weight_step3=`cat $output_dir/final_model_step3/counts_loss_weight.txt`
    echo -e "json_string\t"$step2_bias_json"\nweights\t"$step2_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3/params.txt
    params=$output_dir/final_model_step3/params.txt
   for fold in 0
   do
        min_logcount=4.0
        max_logcount=11.5
        bash $PWD/src/models/chrombpnet_scripts/final_model_step3/train.sh $fold $gpu $model_name $seed $output_dir/final_model_step3 $params  $data_dir/tiledb/db $cell_line $PWD/src/models/chrombpnet_scripts/final_model_step3/bpnet_with_bias.py $min_logcount $max_logcount
        bash $PWD/src/models/chrombpnet_scripts/final_model_step3/predict.sh $fold $gpu $model_name $seed  $output_dir/final_model_step3  $data_dir/tiledb/db $cell_line $chrom_sizes
        bash $PWD/src/models/chrombpnet_scripts/final_model_step3/score.sh $output_dir/final_model_step3 $model_name $fold $cell_line $seed $min_logcount $max_logcount
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3
fi


#bash $PWD/src/models/chrombpnet_scripts/final_model_step3/score.sh $output_dir/final_model_step3 $model_name $fold $cell_line $seed

## UNPLUG MODEL

if [[ -d $output_dir/final_model_step3/unplug ]] ; then
    echo "skipping unplugging"
else
    echo "UNPLUGGING!!!!" 
    mkdir $output_dir/final_model_step3/unplug
    counts_loss_weight_step3=`cat $output_dir/final_model_step3/counts_loss_weight.txt`
    unplug_bias_json=$output_dir/final_model_step3/microglia.0.arch
    unplug_bias_weights=$output_dir/final_model_step3/microglia.0.weights
    echo -e "json_string\t"$unplug_bias_json"\nweights\t"$unplug_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3/unplug/params.txt

    params=$output_dir/final_model_step3/unplug/params.txt

    for fold in 0
    do
        min_logcount=4.0
        max_logcount=11.5
        CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/models/chrombpnet_scripts/unplug/get_model_with_bias_unplugged.py --model_params $params --outf $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 
        bash  $PWD/src/models/chrombpnet_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/final_model_step3/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
        bash  $PWD/src/models/chrombpnet_scripts/unplug/score.sh $output_dir/final_model_step3/unplug $model_name $fold $cell_line $seed $min_logcount $max_logcount
    done
    cp $PWD/$cur_file_name $output_dir/final_model_step3/unplug
fi

min_logcount=4.0
max_logcount=11.5
#bash  $PWD/src/models/chrombpnet_scripts/unplug/score.sh $output_dir/final_model_step3/unplug $model_name $fold $cell_line $seed $min_logcount $max_logcount

### GET INTERPRETATIONS


if [[ -d $output_dir/final_model_step3/unplug/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/final_model_step3/unplug/deepshap
    #bed_file=$data_dir/$cell_line"_idr_split"
    bed_file=data/microglia.idr.optimal.narrowPeak
    for fold in 0
    do
        bash $PWD/src/evaluation/interpret/interpret.sh \
	     $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 \
	     $bed_file \
	     0 \
	     $data_dir/tiledb/db \
	     $chrom_sizes \
	     $output_dir/final_model_step3/unplug/deepshap \
	     $cell_line \
	     $gpu \
	     $fold
	
        #bash $PWD/src/evaluation/interpret/interpret.sh $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3/unplug/deepshap $cell_line $gpu $fold
        #bash $PWD/src/evaluation/interpret/interpret.sh $output_dir/final_model_step3/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3/unplug/deepshap $cell_line $gpu $fold
    done

    #python $PWD/src/evaluation/interpret/combine_shap_pickle.py --source $output_dir/final_model_step3/unplug/deepshap --target $output_dir/final_model_step3/unplug/deepshap --type 20k
    #cp $PWD/$cur_file_name $output_dir/final_model_step3/unplug/deepshap

fi

modisco_sig_dir=modisco #/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/

if [[ -d $modisco_sig_dir/$cell_line/$setting"_new"/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line/$setting"_new"/
    modisco_dir_final=$modisco_sig_dir/$cell_line/$setting"_new"/
    cp  $output_dir/final_model_step3/unplug/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi


#BIAS INTERPRETATIONS


if  [[ -d $output_dir/invivo_bias_model_step1/ ]]  
thenx
    if [[ -d $output_dir/invivo_bias_model_step1/deepshap ]] ; then
        echo "skipping bias interpretations"
    else
        mkdir $output_dir/invivo_bias_model_step1/deepshap
        #bed_file=$data_dir/$cell_line"_idr_split"
	bed_file=data/microglia.idr.optimal.narrowPeak
    
        for fold in 0
        do
            bash $PWD/src/evaluation/interpret/interpret_weight.sh $output_dir/invivo_bias_model_step1/$model_name.$fold $bed_file 0 $data_dir/tiledb/db $chrom_sizes $output_dir/invivo_bias_model_step1/deepshap $cell_line $gpu $fold
            #bash $PWD/src/evaluation/interpret/interpret_weight.sh $output_dir/invivo_bias_model_step1/$model_name.$fold $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/invivo_bias_model_step1/deepshap $cell_line $gpu $fold
            #bash $PWD/src/evaluation/interpret/interpret_weight.sh $output_dir/invivo_bias_model_step1/$model_name.$fold $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/invivo_bias_model_step1/deepshap $cell_line $gpu $fold
        done

        python $PWD/src/evaluation/interpret/combine_shap_pickle.py --source $output_dir/invivo_bias_model_step1/deepshap --target $output_dir/invivo_bias_model_step1/deepshap --type 20k
        cp $PWD/$cur_file_name $output_dir/invivo_bias_model_step1/deepshap
    fi
else
    echo "skipping step1 interpretations - input bias model given"
fi

modisco_bias_dir=modisco_bias #/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/

if [[ -d $modisco_bias_dir/$cell_line/$setting/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_bias_dir/$cell_line/$setting/
    modisco_dir_final=$modisco_bias_dir/$cell_line/$setting/
    cp  $output_dir/invivo_bias_model_step1/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi


## MAKE FOOTPRINTS


#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs tn5_c --model_dir $output_dir/final_model_step3/unplug/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set1 --model_dir $output_dir/final_model_step3/unplug/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set2 --model_dir $output_dir/final_model_step3/unplug/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs tn5 --model_dir $output_dir/invivo_bias_model_step1/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set2 --model_dir $output_dir/invivo_bias_model_step1/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set1 --model_dir $output_dir/invivo_bias_model_step1/
#CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs gm12878_motifs_set2 --model_dir $output_dir/invivo_bias_model_step1/

