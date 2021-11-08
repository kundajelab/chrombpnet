#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=GM12878
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
setting=$data_type"_"$date"_withinvivobias"
cur_file_name="gm12878_atac_with_invivobias_final.sh"
setting=testing

### SIGNAL INPUT
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
is_filtered=True
samtools_flag=None

blacklist_region=$PWD/data/GRch38_unified_blacklist.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$PWD/results/chrombpnet/$data_type/$cell_line/data_new
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

### MODEL PARAMS

gpu=1
seed=1234 
model_name=model 
neg_dir=$main_dir/negatives_data_new
flank_size=1057

bias_n_dil_layers=8
bias_filters=500

n_dil_layers=8
filters=500

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
    bash $PWD/src/utils/make_gc_matched_negatives/run.sh $neg_dir $overlap_peak $ref_fasta $chrom_sizes $flank_size $stride $blacklist_region 
fi


###  BIAS  INPUTS

bias_json=""
bias_weights=""
neg_bed_train=$neg_dir"/bpnet.inputs.train.flank300.0.negatives.bed"
neg_bed_test=$neg_dir"/bpnet.inputs.test.0.negatives.bed"


### CREATE BIGWIGS

if [[ -d $data_dir ]] ; then
    echo "skipping bigwig creation"
else
    mkdir $data_dir
    bash $PWD/src/utils/preprocess.sh $in_bam $data_dir $samtools_flag $is_filtered $data_type $neg_shift
    cp $PWD/$cur_file_name $data_dir
fi


### CREATE FOLDER TILEDB AND RUN TILEDB

if [[ -d $data_dir/tiledb ]] ; then
    echo "skipping tiledb"
else
    mkdir $data_dir/tiledb
    echo -e "dataset\tnegatives_peak\tidr_peak\toverlap_peak\tambig_peak\tcount_bigwig_unstranded_5p\n"$cell_line"\t"$neg_dir/bpnet.inputs.all.negatives.bed"\t"$idr_peak"\t"$overlap_peak"\t"$blacklist_region"\t" > $data_dir/tiledb/inputs.tsv
    echo -e "overlap_peak\tbed_summit_from_last_col\nnegatives_peak\tbed_summit_from_last_col\nidr_peak\tbed_summit_from_last_col\nambig_peak\tbed_no_summit\ncount_bigwig_unstranded_5p\tbigwig" > $data_dir/tiledb/attribs.txt
    bash $PWD/src/utils/db_ingest.sh  $data_dir/tiledb/inputs.tsv $data_dir/tiledb/db $chrom_sizes $data_dir/tiledb/attribs.txt
    cp $PWD/$cur_file_name $data_dir/tiledb
fi

### STEP 1 -  TRAIN BIAS MODEL

fold=0
min_logcount=0.0
max_logcount=10.0
mkdir $output_dir/invivo_bias_model_step1
params=/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/invivo_bias_model_step1/params.txt
#bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/train_testing.sh $fold $gpu $model_name $seed $output_dir/invivo_bias_model_step1 $params  $data_dir/tiledb/db $cell_line $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/bpnet_model.py $neg_bed_train $min_logcount $max_logcount $neg_bed_train $data_dir/shifted.sorted.bam.bpnet.unstranded.bw
CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/training/train.py \
	--genome=$ref_fasta \
	--tdb_array=$data_dir/tiledb/db \
    --tdb_dataset=$cell_line \
    --peaks=$idr_peak \
    --nonpeaks=$neg_dir"/bpnet.inputs.all.negatives.bed" \
    --output_prefix=$output_dir/invivo_bias_model_step1/model.0 \
    --generator=batchgen \
    --fold=0 \
    --epochs=40 \
    --params=$output_dir/params.txt \
    --inputlen=2114 \
    --outputlen=1000 \
    --bigwig=$data_dir/shifted.sorted.bam.bpnet.unstranded.bw \
    --architecture_from_file=$PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/bpnet_model.py \
    --trackables logcount_predictions_loss loss profile_predictions_loss val_logcount_predictions_loss val_loss val_profile_predictions_loss \




