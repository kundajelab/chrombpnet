#intersect GM12878 DNASE idr peaks with idr peaks for TF chipseq in GM12878 (or close cell line) for pioneer TF's that we have bQTL data for + CTCF
#these intersections will be used to label the umap of all the idr peaks 
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/ENCFF113EFE.bed.gz | bedtools sort -i - | cut -f1,2,3 | uniq > peak_intersections/POU2F2.bed
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/ENCFF071ZMW.bed.gz | bedtools sort -i - | cut -f1,2,3 | uniq > peak_intersections/SPI1.bed 
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/ENCFF323QQU.bed.gz | bedtools sort -i - | cut -f1,2,3 | uniq > peak_intersections/STAT1.bed
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/ENCFF873DJD.bed.gz | bedtools sort -i - | cut -f1,2,3 | uniq > peak_intersections/JUND.bed
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/GSM935478_hg19_wgEncodeSydhTfbsGm12878NfkbTnfaIggrabPk.narrowPeak.gz | bedtools sort -i - | cut -f1,2,3 |  uniq > peak_intersections/NFKB.bed
bedtools intersect -a idr.optimal_peak.narrowPeak.gz -b tf_chipseq_peaks/ENCFF960ZGP.bed.gz | bedtools sort -i - | cut -f1,2,3 | uniq > peak_intersections/CTCF.bed


