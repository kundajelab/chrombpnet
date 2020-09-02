import pdb
from kerasAC.generators.tiledb_predict_generator import *
import imp
#for training
#at summits
#gen=TiledbPredictGenerator(ref_fasta="/mnt/data/annotations/by_release/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
#                    batch_size=20,
#                    tdb_indexer="tasks.tsv",
#                    tdb_partition_attribute_for_upsample="idr_peak",
#                    tdb_partition_thresh_for_upsample=2,
#                    tdb_inputs=['seq'],
#                    tdb_input_source_attribute=['seq'],
#                    tdb_input_flank=[6500],
#                    upsample_ratio=1.0,
#                    tdb_outputs=['tasks.tsv'],
#                    tdb_output_source_attribute=['fc_bigwig'],
#                    tdb_output_flank=[1585],
#                    tdb_input_aggregation=[None],
#                    tdb_input_transformation=[None],
#                    tdb_output_aggregation=[None],
#                    tdb_output_transformation=['asinh'],
#                    num_inputs=1,
#                    num_outputs=1,
#                    chroms=['chr1'],
#                    chrom_sizes="/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes",
#                    expand_dims=False)
#
#X,y,coords=gen[300]
#pdb.set_trace() 

gen=TiledbPredictGenerator(ref_fasta="/mnt/data/annotations/by_release/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
                           batch_size=20,
                           tdb_indexer="tasks.tsv",
                           tdb_inputs=['seq'],
                           tdb_input_source_attribute=['seq'],
                           tdb_input_flank=[6500],
                           tdb_partition_attribute_for_upsample=None,
                           tdb_partition_thresh_for_upsample=None,
                           upsample_ratio=None,
                           tdb_outputs=['tasks.tsv'],
                           tdb_output_source_attribute=['fc_bigwig'],
                           tdb_output_flank=[1585],
                           tdb_input_aggregation=[None],
                           tdb_input_transformation=[None],
                           tdb_output_aggregation=[None],
                           tdb_output_transformation=['asinh'],
                           num_inputs=1,
                           num_outputs=1,
                           chroms=['chr1'],
                           chrom_sizes="/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes",
                           tiledb_stride=10000,
                           expand_dims=False)

#X,y,coords=gen[300]
pdb.set_trace() 
