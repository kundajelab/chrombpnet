cell_line=$1
in_bam=$2
merged_overlap=$3

rgt-hint tracks --bc --bigWig --organism=hg38 $in_bam $merged_overlap  --output-prefix=$cell_line/bias_corrected --forward-shift 4 --reverse-shift -4


