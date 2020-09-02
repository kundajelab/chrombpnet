from kerasAC.custom_losses import *
tdb_path="/mnt/data/encode_dnase_tiledb/db/dnase/ENCSR000EOT"
chrom='chr1'
label_attribute='count_bigwig_unstranded_5p'
ambig_attribute='ambig_peak'
upsample_attribute='idr_peak'
tdb_partition_thresh_for_upsample=1

count_loss_weight=get_loss_weights(tdb_path,chrom,label_attribute,ambig_attribute,upsample_attribute,tdb_partition_thresh_for_upsample)
print(count_loss_weight) 
