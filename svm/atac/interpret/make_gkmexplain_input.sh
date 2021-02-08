for cell_line in GM12878 HEPG2 H1ESC IMR90
do
    python make_gkmexplain_input.py -narrowPeak /srv/scratch/annashch/encode_dnase_tiledb/$cell_line.atac.idr.optimal_peak.narrowPeak.gz -ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -n_to_sample 10000 -outf $cell_line.ATAC.svm.10kb.to.interpret.fa -flank 500
done
