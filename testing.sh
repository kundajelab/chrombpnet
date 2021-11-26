#!/bin/bash

cell_line=GM12878
data_type="ATAC_PE"

date=$(date +'%m.%d.%Y')
#setting=$data_type"_"$date"_withinvivobias"
setting=ATAC_PE_11.15.2021_withinvivobias
cur_file_name="testing.sh" # script used to run for checkpoint

## SIGNAL INPUT
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/GM12878.atac.filt.merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz

## REFERENCE INPUT
blacklist_region=/mnt/data/annotations/blacklist/GRch38/GRch38_unified_blacklist.bed
chrom_sizes=/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes
reference_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

## OUTPUT DIRECTORY FOR STORAGE - Make sure main_dir already exists
main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting
negatives_dir=$main_dir/negatives_data
bigwig_dir=$PWD/results/chrombpnet/$data_type/$cell_line/data

## MODEL PARAMS 
gpu=3
seed=1234 
inputlen=2114
outputlen=1000
fold=0

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

## MAKE NEGATIVES THAT MATCH OVERLAP PEAK SET 

### To build the exlusion bed we will use the overlap peak as they are less conservative compared to idr peak set. This will give as  strict negative set.
### We will also slop both our blacklist and overlap peak set with the flank size and merge them. We will use this combined set as exlusion bed.

if [[ -d $negatives_dir ]] ; then
    echo "negatives director already exists"
else
    mkdir $negatives_dir
    # create an exlusion bed list which is very strict
    flank_size=$(( inputlen/2 ))
    bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $negatives_dir/blacklist_slop1057.bed
    zcat $overlap_peak | awk -v OFS="\t" '{print $1,$2+$10-$flank_size,$2+$10+$flank_size}' > $negatives_dir/peaks_slop1057.bed
    cat $negatives_dir/blacklist_slop1057.bed $negatives_dir/peaks_slop1057.bed | bedtools sort | bedtools merge -i stdin > $negatives_dir/exclude.bed
    exclude_bed=$negatives_dir/exclude.bed

    # create gc-matched negatives
    genomewide_gc=/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_$inputlen.bed
    bash $PWD/src/helpers/make_gc_matched_negatives/run.sh $overlap_peak $exclude_bed $inputlen $negatives_dir $reference_fasta $genomewide_gc
    awk -v OFS="\t" '{print $1, $2, $3, ".",  ".", ".", ".", ".", ".", "500"}' $negatives_dir/negatives.bed > $negatives_dir/negatives_with_summit.bed
    rm $negatives_dir/blacklist_slop1057.bed
    rm $negatives_dir/peaks_slop1057.bed
    # checkpoint cuurent script 
    cp $PWD/$cur_file_name $negatives_dir   
fi

## MAKE BIGWIGS FROM BAM FILES

if [[ -d $bigwig_dir ]] ; then
    echo "skipping bigwig creation"
else
    mkdir $bigwig_dir
    bash $PWD/src/helpers/preprocessing/bam_to_bigwig.sh $in_bam $bigwig_dir $data_type $chrom_sizes 
    # generate pwm matrix from bigwig
    python $PWD/src/helpers/preprocessing/analysis/build_pwm_from_bigwig.py -i $bigwig_dir/shifted.sorted.bam.chrombpnet.unstranded.bw  -g $reference_fasta -o $bigwig_dir/bias_pwm.png -c "chr20" -cz $chrom_sizes 
    # NOTE: open the bias_pwm.png created in the $bigwig_dir to check if the bias pwm looks correct
    # checkpoint current script 
    cp $PWD/$cur_file_name $bigwig_dir 
fi

## GENERATE CHROMBPNET HYPER-PARAMETERS

# edit the parameter files generated from this if needed
# make sure the bias model is trained on atleast 15K regions - adjuts threshold accordingly - after training bias model make sure there are no TF motifs captured by the bias model - do modisco on the peaks

#python $PWD/src/helpers/hyperparameters/find_base_hyperparams.py -b $bigwig_dir/shifted.sorted.bam.chrombpnet.unstranded.bw -p $overlap_peak -n $negatives_dir/negatives_with_summit.bed \
#        -g $reference_fasta -sr 0.1 -t 0.8  \
#        -fl $fold \
#        -ol 1000 \
#        -o $output_dir 

## STEP 1 - TRAIN and TEST BIAS MODEL

mkdir $output_dir/invivo_bias_model_step1

# CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/training/train.py \
#         --genome=$reference_fasta \
#         --peaks=$overlap_peak \
#         --nonpeaks=$negatives_dir"/negatives_with_summit.bed" \
#         --output_prefix=$output_dir/invivo_bias_model_step1/model.0 \
#         --generator=batchgen \
#         --fold=$fold \
#        --epochs=40 \
#        --params=$output_dir/params.txt \
#        --inputlen=2114 \
#        --outputlen=1000 \
#        --max-jitter=10 \
#        --batch-size=64 \
#        --bigwig=$bigwig_dir/shifted.sorted.bam.chrombpnet.unstranded.bw \
#        --architecture_from_file=$PWD/src/training/models/bpnet_model.py \
#        --trackables logcount_predictions_loss loss logits_profile_predictions_loss val_logcount_predictions_loss val_loss val_logits_profile_predictions_loss 

# CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/training/predict.py \
#         --genome=$reference_fasta \
#         --peaks=$overlap_peak \
#         --nonpeaks=$negatives_dir"/negatives_with_summit.bed" \
#         --output_prefix=$output_dir/invivo_bias_model_step1/ \
#         --generator=batchgen \
#         --fold=$fold \
#         --inputlen=2114 \
#         --outputlen=1000 \
#         --batch-size=64 \
#         --model_h5=$output_dir/invivo_bias_model_step1/model.0.h5 \
#         --bigwig=$bigwig_dir/shifted.sorted.bam.chrombpnet.unstranded.bw \

## STEP 2 - TRAIN SCALED BIAS MODEL


## STEP 3 - TRAIN and TEST BIAS MODEL WITH CHROMBPNET


## INTERPRETING

if [[ -d $main_dir/"subsample_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $main_dir/"subsample_idr_split" 
    flank_size=$(( inputlen/2 ))
    bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > temp.txt
    bedtools intersect -v -a $idr_peak -b temp.txt | shuf  > $main_dir/"subsample_idr_split/temp.txt"
    split -l 10000 $main_dir/"subsample_idr_split/temp.txt" $main_dir/"subsample_idr_split/x"
    cat $main_dir/"subsample_idr_split/xaa" $main_dir/"subsample_idr_split/xab" $main_dir/"subsample_idr_split/xac" > $main_dir/"subsample_idr_split/30K.subsample.idr.bed"
    rm  $main_dir/"subsample_idr_split/temp.txt"
    rm temp.txt
fi


CUDA_VISIBLE_DEVICES=$gpu python $PWD/src/evaluation/interpret/interpret.py \
         --genome=$reference_fasta \
         --regions=$main_dir/"subsample_idr_split/30K.subsample.idr.bed"\
         --output_prefix=$main_dir/"subsample_idr_split/30K.subsample.idr.bed" \
         --model_h5=$output_dir/invivo_bias_model_step1/model.0.h5 \

## MARGINAL FOOTPRINTING 


## INTERPRET MODELS ON RANDOM SUBSET OF REGIONS - USED FOR MODISCO














