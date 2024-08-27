modeldir=$1
celll=$2
foldn=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits_format/$3.json
tag=$4
fold=$5

python get_gc_matched_negatives_test.py \
	-c /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$tag/$celll/negatives_data$fold/candidate.negatives.bed \
	-f /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$tag/$celll/negatives_data$fold/foreground.gc.bed \
	-o $modeldir/train_test_regions_may_7_2024/negatives \
	-fl $foldn		

awk -v OFS="\t" '{print $1, $2, $3, ".",  ".", ".", ".", ".", ".", "1057"}' $modeldir/train_test_regions_may_7_2024/negatives.bed > $modeldir/train_test_regions_may_7_2024/negatives_with_summit.bed

gzip -c  $modeldir/train_test_regions_may_7_2024/negatives_with_summit.bed >  $modeldir/train_test_regions_may_7_2024/nonpeaks.testset.bed.gz

