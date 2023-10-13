bw1=results/chrombpnet/ATAC_PE/K562/data/K562_unstranded.bw
bw2=results/chrombpnet/ATAC_PE/K562/data/K562_unstranded.bw
output=temps.bw
wiggletools mean $bw1 $bw2 | wigToBigWig stdin /mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes $output

