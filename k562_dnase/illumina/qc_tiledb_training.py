import pdb
from kerasAC.generators.tiledb_generator import *
import imp
#for training
gen=TiledbGenerator(ref_fasta="/mnt/data/annotations/by_release/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
                    batch_size=20,
                    tdb_indexer="tasks.tsv",
                    tdb_partition_attribute_for_upsample="idr_peak",
                    tdb_partition_thresh_for_upsample=1,
                    tdb_inputs=['seq'],
                    tdb_input_source_attribute=['seq'],
                    tdb_input_flank=[6500],
                    upsample_ratio=0.3,
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
                    shuffle_epoch_start=True,
                    shuffle_epoch_end=True,
                    pseudocount=0,
                    add_revcomp=True,
                    expand_dims=False,
                    return_coords=True)



#load the model 
architecture_module=imp.load_source('','test_model.py')
model=architecture_module.getModelGivenModelOptionsAndWeightInits()
for i in range(20,200):
    X,y,coords=gen[i]
    model.train_on_batch(X,y)

X,y,coords=gen[300]
preds=model.predict(X)
print(preds)
pdb.set_trace() 
