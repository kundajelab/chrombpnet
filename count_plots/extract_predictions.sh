#K562 DNASE 
#corrected vs uncorrected 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/seed.1234.cw.250.filters.500.naive.range.4.6.to.11.5/predictions.k562.withdups.1234seed.250counts.$fold.predictions.1 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/with_bias_bpnet_tobias/k562.dnase.with.bpnet.tobias.bias.$fold.predictions \
	   -formats legacy new \
	   -labels K562_DNASE_Uncorrected.$fold K562_DNASE_Corrected_BPNET_TOBIAS.$fold \
	   -outf K562_DNASE_Uncorrected_vs_Corrected.$fold
done
#corrected vs unplugged 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.$fold.predictions /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/with_bias_bpnet_tobias/k562.dnase.with.bpnet.tobias.bias.$fold.predictions \
	   -formats new new \
	   -labels K562_DNASE_Unplugged.$fold K562_DNASE_Corrected_BPNET_TOBIAS.$fold \
	   -outf K562_DNASE_Unplugged_vs_Corrected.$fold

done

#HEPG2 DNASE
#corrected vs uncorrected 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/seed.1234.cw.101.filters.500.naive.range.4.6.to.11.5/seed.1234.cs.101.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions.1 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/with_bias_bpnet_tobias/hepg2.dnase.with.bpnet.tobias.bias.$fold.predictions \
	   -formats legacy new \
	   -labels HEPG2_DNASE_Uncorrected.$fold HEPG2_DNASE_Corrected_BPNET_TOBIAS.$fold \
	   -outf HEPG2_DNASE_Uncorrected_vs_Corrected.$fold
done
#corrected vs unplugged 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/with_bias_unplugged/hepg2.dnase.with.bias.unplugged.$fold.predictions /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/with_bias_bpnet_tobias/hepg2.dnase.with.bpnet.tobias.bias.$fold.predictions \
	   -formats new new \
	   -labels HEPG2_DNASE_Unplugged.$fold HEPG2_DNASE_Corrected_BPNET_TOBIAS.$fold \
	   -outf HEPG2_DNASE_Unplugged_vs_Corrected.$fold
done

#HEPG2 ATAC
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/hepg2_atac/baseline/hepg2.atac.seed.1234.cs.50.filters.500.naive.range.4.6.to.11.5.$fold.predictions.1 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_tobias/hepg2.atac.with.tobias.bias.$fold.predictions \
	   -formats legacy new \
	   -labels HEPG2_ATAC_Uncorrected.$fold HEPG2_ATAC_Corrected_TOBIAS.$fold \
	   -outf HEPG2_ATAC_Uncorrected_vs_Corrected.$fold
done
#corrected vs unplugged 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_unplugged/hepg2.atac.with.bias.unplugged.$fold.predictions /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_tobias/hepg2.atac.with.tobias.bias.$fold.predictions \
	   -formats new new \
	   -labels HEPG2_ATAC_Unplugged.$fold HEPG2_ATAC_Corrected_TOBIAS.$fold \
	   -outf HEPG2_ATAC_Unplugged_vs_Corrected.$fold
done

#GM12878 ATAC
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/gm12878_atac/baseline/gm12878.atac.seed.1234.cs.55.filters.500.naive.range.4.6.to.11.5.$fold.predictions.1 /srv/scratch/annashch/deeplearning/profile/gm12878_atac/with_bias_tobias/gm12878.atac.with.tobias.bias.$fold.predictions \
	   -formats legacy new \
	   -labels GM12878_ATAC_Uncorrected.$fold GM12878_ATAC_Corrected_TOBIAS.$fold \
	   -outf GM12878_ATAC_Uncorrected_vs_Corrected.$fold
done
#corrected vs unplugged 
for fold in `seq 0 4`
do
    python extract_predictions.py -predictions_files /srv/scratch/annashch/deeplearning/profile/gm12878_atac/with_bias_unplugged/gm12878.atac.with.bias.unplugged.$fold.predictions /srv/scratch/annashch/deeplearning/profile/gm12878_atac/with_bias_tobias/gm12878.atac.with.tobias.bias.$fold.predictions \
	   -formats new new \
	   -labels GM12878_ATAC_Unplugged.$fold GM12878_ATAC_Corrected_TOBIAS.$fold \
	   -outf GM12878_ATAC_Unplugged_vs_Corrected.$fold
done
