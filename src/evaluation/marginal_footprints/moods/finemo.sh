peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed
fasta=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
bigwigs=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/interpret/merged.GM12878.counts.bw
output=finemo_output/fold0
modisco=/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/SIGNAL/modisco_crop_500/modiscolite_counts.h5
outputp=$output/gm.atac.572M.extracts
outputd=$output

#finemo extract-regions-bw  -p $peaks -f $fasta -b $bigwigs -o $outputp -w 1000

outoutp=$output/gm.atac.572M.extracts.npz

CUDA_VISIBLE_DEVICES=3 finemo call-hits -r $outoutp -m $modisco -o $outputd -p $peaks

finemo report -r $outoutp -H $outputd/hits_unique.tsv -p $peaks -m $modisco -o $outputd
