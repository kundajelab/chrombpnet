#awk 'BEGIN {OFS="\t"} {print $2, $3,$4, $1, $5, $6, $7, $8, $9}' chr1_fimo.tsv | grep -v 'fimo' > formatted_input.bed

bedtools slop -i ../merged_output.bed -g /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes -b 100 > inpeaks_hits.bed

bedtools intersect -v -a formatted_input.bed -b inpeaks_hits.bed -wa > unbound_hits.tsv 


