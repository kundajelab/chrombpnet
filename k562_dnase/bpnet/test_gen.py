import kerasAC 
from kerasAC.generators.tiledb_generator import *
from kerasAC.tiledb_config import *
import tiledb 
import pdb

gen=TiledbGenerator(ref_fasta="/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
                    batch_size=50,
                    tdb_array="/srv/scratch/annashch/encode_dnase_tiledb/db_singlethread/dnase",
                    tdb_partition_attribute_for_upsample="idr_peak",
                    tdb_partition_thresh_for_upsample=1,
                    tdb_input_flank=[673],
                    tdb_input_source_attribute=["seq"],
                    tdb_input_aggregation=[None],
                    tdb_input_transformation=[None],
                    tdb_output_source_attribute=["count_bigwig_unstranded_5p","count_bigwig_unstranded_5p"],
                    tdb_output_flank=[500,500],
                    tdb_output_aggregation=[None,"sum"],
                    tdb_output_transformation=[None,"asinh"],
                    num_inputs=1,
                    num_outputs=2,
                    chroms=['chr2','chr3','chr4','chr5','chr6','chr7','chr9','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'],
                    upsample_ratio=None,
                    tasks=['ENCSR000EOT'],
                    tdb_ambig_attribute='ambig_peak',
                    shuffle_epoch_start=True,
                    shuffle_epoch_end=True,
                    return_coords=True,
                    tdb_config=get_default_config(),
                    tdb_ctx=tiledb.Ctx(config=get_default_config()))
pdb.set_trace()
