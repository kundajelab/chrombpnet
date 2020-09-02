for f in `cat atac.tier1.overlap.peak.files.txt`
do
    echo $f
    zcat $f >> atac.tier1.overlap.concat.bed  
done
bedtools sort -i atac.tier1.overlap.concat.bed > atac.tier1.overlap.concat.sorted.bed
bedtools merge -i atac.tier1.overlap.concat.sorted.bed > atac.tier1.overlap.merged.bed


    
