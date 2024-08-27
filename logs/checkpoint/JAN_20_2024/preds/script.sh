modelpath=$1
regions=$2
genom=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
chroms=/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes
outpref=$3
gpu=$4

CUDA_VISIBLE_DEVICES=$gpu  python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/make_bigwigs/predict_to_bigwig_new.py \
	-cm $modelpath/new_model_formats/chrombpnet/ \
	-cmb $modelpath/new_model_formats/chrombpnet_wo_bias/ \
	-r $regions \
	-g $genom \
	-c $chroms \
	-o $outpref 

