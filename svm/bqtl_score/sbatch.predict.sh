#!/bin/bash
atac_model=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/atac/GM12878/models/model.GM12878.0.txt
dnase_model=/oak/stanford/groups/akundaje/projects/chrombpnet/svm/dnase/GM12878/models/GM12878.0.model.txt
for tf in nfkB junD oct11 oct1 pu1 stat1 
do
    for allele in ref alt
    do
	sbatch -J $tf.$allele -o logs/$tf.$allele.o -e logs/$tf.$allele.e --mem=30000 --time=1440 --mincpus=16 -p akundaje,euan,normal predict.sh /oak/stanford/groups/akundaje/projects/chrombpnet/svm/bQTL/gm12878_gkm_preds_bqtl_$tf.$allele.fa $atac_model /oak/stanford/groups/akundaje/projects/chrombpnet/svm/bQTL/preds_$allele/preds.gm12878.0.$tf.$allele
	sbatch -J $tf.$allele -o logs/$tf.$allele.o -e logs/$tf.$allele.e --mem=30000 --time=1440 --mincpus=16 -p akundaje,euan,normal predict.sh /oak/stanford/groups/akundaje/projects/chrombpnet/svm/bQTL/gm12878_gkm_preds_bqtl_$tf.$allele.fa $atac_model /oak/stanford/groups/akundaje/projects/chrombpnet/svm/bQTL/preds_atac_$allele/preds.gm12878.0.$tf.$allele
    done
done
