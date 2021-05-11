for f in /mnt/lab_data3/anusri/spi1_deepshaps/all_deepshap/H3K4me3.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/all_deepshap/gm12878.dnase.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/all_deepshap/H3K27Ac.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/all_deepshap/H3K4me1.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k27ac.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.dnase.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k4me1.strand1.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k4me1.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k4me3.strand0.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k27ac.strand1.fold0.deepSHAP /mnt/lab_data3/anusri/spi1_deepshaps/tf_deepshap/SPI1.h3k4me3.strand1.fold0.deepSHAP 
do
python deepshap_to_bigwig.py --path_to_pickle $f \
       --output_dir /srv/scratch/annashch/chrombpnet/make_tracks_from_anusri_deepshap  &
done

