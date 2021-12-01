reference_fasta=$1
chrom_sizes=$2
blacklist_region=$3
overlap_peak=$4
inputlen=$5
output_dir=$6
genomewide_gc=$7

# To build the exlusion bed we will use the overlap peak as they are less conservative compared to idr peak set. This will give as  strict negative set.
# We will also slop both our blacklist and overlap peak set with the inputlen//2 and merge them. We will use this combined set as exlusion bed.
# we slop to make sure that the regions we choose have no intersection (not even 1bp) with the overlap peaks/blacklist regions
flank_size=$(( inputlen/2 ))
bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > $output_dir/blacklist_slop1057.bed
zcat $overlap_peak | awk -v OFS="\t" '{print $1,$2+$10-$flank_size,$2+$10+$flank_size}' > $output_dir/peaks_slop1057.bed
cat $output_dir/blacklist_slop1057.bed $output_dir/peaks_slop1057.bed | bedtools sort | bedtools merge -i stdin > $output_dir/exclude.bed
rm $output_dir/blacklist_slop1057.bed
rm $output_dir/peaks_slop1057.bed
exclude_bed=$output_dir/exclude.bed

# create regions of size inputlen that do not fall in exclude bed and that gc-match with the given overlap peaks
# create as many regions as there are in overlap peaks bed file
bash $PWD/src/helpers/make_gc_matched_negatives/run.sh $overlap_peak $exclude_bed $inputlen $output_dir $reference_fasta $genomewide_gc
# make a dummy summit file - rest of the pipeline uses this format
awk -v OFS="\t" '{print $1, $2, $3, ".",  ".", ".", ".", ".", ".", "500"}' $output_dir/negatives.bed > $output_dir/negatives_with_summit.bed

