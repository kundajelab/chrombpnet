#python make_gkmexplain_input.py -narrowPeak /srv/scratch/annashch/encode_dnase_tiledb/K562.dnase.idr.optimal_peak.narrowPeak.gz -ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -n_to_sample 10000 -outf K562.DNASE.svm.10kb.to.interpret.fa -flank 500

#python make_gkmexplain_input.py -narrowPeak /srv/scratch/annashch/encode_dnase_tiledb/HEPG2.dnase.idr.optimal_peak.narrowPeak.gz -ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -n_to_sample 10000 -outf HEPG2.DNASE.svm.10kb.to.interpret.fa -flank 500


for cell_line in  GM12878 H1ESC IMR90
do
    python make_gkmexplain_input.py -narrowPeak /srv/scratch/annashch/encode_dnase_tiledb/$cell_line.dnase.idr.optimal_peak.narrowPeak.gz -ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -n_to_sample 10000 -outf $cell_line.DNASE.svm.10kb.to.interpret.fa -flank 500
done
