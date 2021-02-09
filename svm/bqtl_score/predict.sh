#!/bin/bash
gkmpredict -v 2 -T 16 $1 $2 $3



#!/bin/bash
#get gkm predictions across all folds for the saturation mutagenesis paper regions
atac_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/atac
dnase_prefix=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/dnase
#HEPG2
for fold in `seq 0 9`
do
    #get predictions
    gkmpredict -v 2 -T 16 regions.hepg2.txt $atac_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.atac.preds.$fold
    gkmpredict -v 2 -T 16 regions.hepg2.txt $dnase_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.dnase.preds.$fol\
d
    #get interpretations
    gkmexplain regions.hepg2.txt $atac_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.atac.gkmexplain.$fold
    gkmexplain regions.hepg2.txt $dnase_prefix/HEPG2/models/HEPG2.$fold.model.txt hepg2.dnase.gkmexplain.$fold
done
