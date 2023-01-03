mode=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/predictions_at_loci_flank_150/jaspar/
mkdir $mode

fasta_in=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/predictions_at_loci_flank_150/sequence.fasta
#for pfm in /oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/THP1/ATAC_PE_03.24.2022_withgm12878bias/SIGNAL/ppms/counts/*
for pfm in /oak/stanford/groups/akundaje/soumyak/motifs/pfms/*
    do
        base=$(basename -- "$pfm")
        tf="${base%.*}"
        tf="${tf%.*}"
        tf="${tf%.*}"
        echo $tf
        moods-dna.py -m $pfm -s $fasta_in -p 0.0001  --batch > $mode/$tf.hits.csv &
   done


wait

#cat $mode/ref/metacluster*.hits.csv > $mode/ref/all_metaclusters.all_patterns.hits.csv
cat $mode/*.hits.csv > $mode/all_patterns.hits.csv




