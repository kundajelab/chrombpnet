input_name=$1
chrom_sizes=$2

python $PWD/main_scripts/make_gc_matched_negatives/make_bed_files.py --input_name $input_name

all_bed=$input_name".all.negatives.bed"
grep -w "chr1" $all_bed > $input_name".test.0.negatives.bed"
grep -w -v "chr1" $all_bed > $input_name".train.0.negatives.bed"
bedtools slop -i  $input_name".train.0.negatives.bed" -g $chrom_sizes -b 300 > $input_name".train.flank300.0.negatives.bed"
