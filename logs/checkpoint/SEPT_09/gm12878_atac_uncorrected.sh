#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=GM12878
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
setting=4_$neg_shift"_shifted_"$data_type"_uncorrected_"$date
cur_file_name="gm12878_atac_uncorrected.sh"

### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak
#.gz file?
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
is_filtered=True
samtools_flag=None

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



### STEP 1 -  TRAIN BIAS MODEL


if [[ -d $output_dir/uncorrected_model ]] ; then
    echo "skipping step 1 - directory present "
else

    if test -z "$bias_json" 
    then
    	mkdir $output_dir/uncorrected_model

        bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/uncorrected_model/counts_loss_weight.txt
        counts_loss_weight_step1=`cat $output_dir/uncorrected_model/counts_loss_weight.txt`
 
   	echo -e "counts_loss_weight\t"$counts_loss_weight_step1"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/uncorrected_model/params.txt
    	params=$output_dir/uncorrected_model/params.txt
    	for fold in 0
    	do
            ./hint_atac_scripts/main_scripts/model/train.sh $fold $gpu $model_name $seed $output_dir/uncorrected_model $params  $data_dir/tiledb/db $cell_line $PWD/hint_atac_scripts/main_scripts/model/profile_bpnet_dnase_with_bias.py
            ./hint_atac_scripts/main_scripts/model/predict.sh $fold $gpu $model_name $seed  $output_dir/uncorrected_model  $data_dir/tiledb/db $cell_line $chrom_sizes
            ./hint_atac_scripts/main_scripts/model/score.sh $output_dir/uncorrected_model $model_name $fold $cell_line $seed
    	done
    	cp $PWD/$cur_file_name $output_dir/uncorrected_model
    else
        echo "skipping step1 - input bias model given"
    fi

fi




##BIAS INTERPRETATIONS


if  [[ -d $output_dir/uncorrected_model/ ]]  
then
    if [[ -d $output_dir/uncorrected_model/deepshap ]] ; then
        echo "skipping bias interpretations"
    else
        mkdir $output_dir/uncorrected_model/deepshap
        bed_file=$data_dir/$cell_line"_idr_split"

        for fold in 0
        do
            ./main_scripts/interpret/interpret_weight.sh $output_dir/uncorrected_model/$model_name.$fold $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/uncorrected_model/deepshap $cell_line $gpu $fold
            ./main_scripts/interpret/interpret_weight.sh $output_dir/uncorrected_model/$model_name.$fold $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/uncorrected_model/deepshap $cell_line $gpu $fold
            ./main_scripts/interpret/interpret_weight.sh $output_dir/uncorrected_model/$model_name.$fold $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/uncorrected_model/deepshap $cell_line $gpu $fold
        done

        python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/uncorrected_model/deepshap --target $output_dir/uncorrected_model/deepshap --type 20k
        cp $PWD/$cur_file_name $output_dir/uncorrected_model/deepshap
    fi
else
    echo "skipping step1 interpretations - input bias model given"
fi



modisco_bias_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/
if [[ -d $modisco_bias_dir/$cell_line ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_bias_dir/$cell_line
fi

if [[ -d $modisco_bias_dir/$cell_line/$setting/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_bias_dir/$cell_line/$setting/
    modisco_dir_final=$modisco_bias_dir/$cell_line/$setting/
    cp  $cell_line/$setting/uncorrected_model/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi

### RUN MODISCO



