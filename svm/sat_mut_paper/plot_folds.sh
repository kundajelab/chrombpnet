#!/bin/bash 
python plot_folds.py --bed regions.hepg2.txt \
       --scores_prefix hepg2.dnase.gkmexplain \
       --outf_prefix hepg2.dnase.gkmexplain.plot
python plot_folds.py --bed regions.hepg2.txt \
       --scores_prefix hepg2.atac.gkmexplain \
       --outf_prefix hepg2.atac.gkmexplain.plot
python plot_folds.py --bed regions.k562.txt \
       --scores_prefix k562.dnase.gkmexplain \
       --outf_prefix k562.dnase.gkmexplain.plot

