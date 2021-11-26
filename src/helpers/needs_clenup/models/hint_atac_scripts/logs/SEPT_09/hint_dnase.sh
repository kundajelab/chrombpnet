
cell_line=$1
in_bam=$2
merged_overlap=$3

rgt-hint tracks --bc --bigWig --organism=hg38 $in_bam $merged_overlap  --bias-table /users/anusri/rgtdata/fp_hmms/double_hit_bias_table_F.txt,/users/anusri/rgtdata/fp_hmms/double_hit_bias_table_R.txt --output-prefix=$PWD/hint_atac_scripts/$cell_line/bias_corrected_dnase --forward-shift 0 --reverse-shift +1


