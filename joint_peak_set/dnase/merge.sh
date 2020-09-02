for f in `cat dnase.tier1.overlap.peak.files.txt`
do
    echo $f
    zcat $f >> dnase.tier1.overlap.concat.bed  
done
bedtools sort -i dnase.tier1.overlap.concat.bed > dnase.tier1.overlap.concat.sorted.bed
bedtools merge -i dnase.tier1.overlap.concat.sorted.bed > dnase.tier1.overlap.merged.bed


    
