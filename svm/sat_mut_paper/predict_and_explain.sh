#!/bin/bash
#get gkm predictions across all folds for the saturation mutagenesis paper regions
atac_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/atac
dnase_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/dnase
#HEPG2 
for fold in `seq 0 9`
do
    #get predictions 
    gkmpredict -v 2 -T 16 hepg2.fasta $atac_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.atac.preds.$fold
    gkmpredict -v 2 -T 16 hepg2.fasta $dnase_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.dnase.preds.$fold
    #get interpretations
    gkmexplain hepg2.fasta $atac_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.atac.gkmexplain.$fold
    gkmexplain hepg2.fasta $dnase_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.dnase.gkmexplain.$fold    
done

#K562
for fold in `seq 0 9`
do
    #get predictions 
    gkmpredict -v 2 -T 16 k562.fasta $dnase_prefix/K562/models/K562.$fold.model.txt k562.dnase.preds.$fold
    #get interpretations
    gkmexplain k562.fasta $dnase_prefix/K562/models/K562.$fold.model.txt k562.dnase.gkmexplain.$fold    
done
