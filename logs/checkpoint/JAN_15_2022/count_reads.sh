in_bam=$1
samtools view -b -@50 -F780 -q30 $in_bam | samtools view -c 
