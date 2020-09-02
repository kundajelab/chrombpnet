import argparse
from kerasAC.architectures.profile_bpnet_dnase import *
from kerasAC.generators.tiledb_generator import *
import pdb

parser=argparse.ArgumentParser(description="view model arch")
parser.add_argument("--seed",type=int,default=1234)
parser.add_argument("--init_weights",default=None)
parser.add_argument("--tdb_input_flank",nargs="+",default=[673])
parser.add_argument("--tdb_output_flank",nargs="+",default=[500])
parser.add_argument("--num_tasks",type=int,default=1)

args=parser.parse_args()
model=getModelGivenModelOptionsAndWeightInits(args)
print("got model") 
gen=TiledbGenerator("/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
                    batch_size=20,
                    tdb_indexer="tasks.tsv",
                    tdb_partition_attribute_for_upsample="idr_peak",
                    tdb_partition_thresh_for_upsample=1,
                    tdb_inputs=['seq'],
                    tdb_input_source_attribute=["seq"],
                    tdb_input_aggregation=[None],
                    tdb_input_transformation=[None],
                    tdb_input_flank=[673],
                    tdb_outputs=['tasks.tsv','tasks.tsv'],
                    tdb_output_source_attribute=['count_bigwig_unstranded_5p','count_bigwig_unstranded_5p'],
                    tdb_output_flank=[500,500],
                    tdb_output_aggregation=[None,'sum'],
                    tdb_output_transformation=[None,'asinh'],
                    num_inputs=1,
                    num_outputs=2,
                    chroms=['chr21'],
                    chrom_sizes="/mnt/data/annotations/by_release/hg38/hg38.chrom.sizes",
                    upsample_ratio=1.0)
print("got generator")
X,y=gen[0]
model.fit(x=X,y=y)
preds=model.predict(X)
pdb.set_trace()

