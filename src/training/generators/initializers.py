import generators.tiledb_generator as tiledb_generator
import generators.batchgen_generator as batchgen_generator
import utils.batchgen_generator_utils as data_utils
import pandas as pd
import splits


# batchgen initializers

def initialize_generators(args, parameters, mode, generator_name):

    bed_regions, chroms=splits.get_bed_regions_for_fold_split(bed_regions, args.fold, mode) 

    # need generator to crop and revcomp aug training examples, but not for 
    # validation. Also applies bias model to cropped, rev comp-ed seqs.
    train_generator = batchgen_generator.ChromBPNetBatchGenerator(bed_regions args.negative_sampling_ratio, 
                                        input_len, output_len, args.batch_size)

    return train_generator

## tiledb initializers 
def initialize_generators_tiledb(args, parameters, mode):

    assert(args.inputlen%2==0)
    assert(args.outputlen%2==0)
    
    # get only those peak/non peak regions corresponding to train/valid/test set
    if args.bed_regions is not None:
            bed_regions=pd.read_csv(args.bed_regions,header=None,sep='\t')
            bed_regions, chroms=splits.get_bed_regions_for_fold_split(bed_regions, args.fold, mode)
    else: 
            bed_regions=None

    if args.nonpeaks is not None:
            nonpeak_regions=pd.read_csv(args.nonpeaks,header=None,sep='\t')
            nonpeak_regions, chroms=splits.get_bed_regions_for_fold_split(nonpeak_regions, args.fold, mode) 
    else:
            nonpeak_regions=None

    if mode=="test":
            # make sure there is no jitter at test time - we are testing only at summits
            assert(args.max_jitter==0)
  
    generator=tiledb_generator.TiledbGenerator(chroms=chroms,
                                    ref_fasta=args.genome,
                                    batch_size=args.batch_size,
                                    tdb_array=args.tdb_array,
                                    tdb_input_source_attribute=args.tdb_input_source_attribute,
                                    tdb_input_flank=args.input_len//2,
                                    pseudocount=args.pseudocount,
                                    tdb_output_source_attribute=args.tdb_output_source_attribute,
                                    tdb_output_flank=args.output_len//2,
                                    tdb_output_min=["None", parameters["count_min_threshold"]],
                                    tdb_output_max=["None", parameters["count_max_threshold"]],
                                    tdb_output_aggregation=["None", "sum"],
                                    tdb_output_transformation=["None", "log"],
                                    tdb_input_datasets=args.tdb_input_datasets,
                                    tdb_output_datasets=args.tdb_output_datasets,
                                    num_inputs=1,
                                    num_outputs=2,
                                    add_revcomp=True,
                                    bed_regions=bed_regions,
                                    nonpeak_regions=nonpeak_regions,
                                    max_jitter=args.max_jitter,
                                    negative_sampling_ratio=args.negative_sampling_ratio
                                    )
    
    print("data generator is ready!")

    return generator
  
