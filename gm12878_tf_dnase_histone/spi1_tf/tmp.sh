#!/bin/bash
module load ucsc_tools

#get the data
#wget https://www.encodeproject.org/files/ENCFF724WMD/@@download/ENCFF724WMD.bam -O rep1.bam
#wget https://www.encodeproject.org/files/ENCFF482TVZ/@@download/ENCFF482TVZ.bam -O rep2.bam
#wget https://www.encodeproject.org/files/ENCFF857FLV/@@download/ENCFF857FLV.bam -O control.bam


#merge
#samtools merge merged.bam rep1.bam rep2.bam

#generate bedgraph and bw file for positive strand
#bedtools genomecov -5 -bg -strand + -g hg38.chrom.sizes -ibam merged.bam | sort -k1,1 -k2,2n > pos_strand.bedGraph
bedGraphToBigWig pos_strand.bedGraph hg38.chrom.sizes pos_strand.bw

#generate bedgraph file for negative strand
#bedtools genomecov -5 -bg -strand - -g hg38.chrom.sizes -ibam merged.bam | sort -k1,1 -k2,2n > neg_strand.bedGraph
bedGraphToBigWig neg_strand.bedGraph hg38.chrom.sizes neg_strand.bw

#generate bedgraph and bw file for control
#bedtools genomecov -5 -bg -strand + -g hg38.chrom.sizes -ibam control.bam | sort -k1,1 -k2,2n > control_pos_strand.bedGraph
bedGraphToBigWig control_pos_strand.bedGraph hg38.chrom.sizes control_pos_strand.bw
#bedtools genomecov -5 -bg -strand - -g hg38.chrom.sizes -ibam control.bam | sort -k1,1 -k2,2n > control_neg_strand.bedGraph
bedGraphToBigWig control_neg_strand.bedGraph hg38.chrom.sizes control_neg_strand.bw
