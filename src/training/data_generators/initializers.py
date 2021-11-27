import generators.tiledb_generator as tiledb_generator
import generators.batchgen_generator as batchgen_generator
import utils.batchgen_generator_utils as data_utils
import pandas as pd
import splits

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def fetch_values_based_on_mode(mode, args, parameters, nonpeak_regions):
    if mode=="train":
        cts_sum_min_thresh=parameters["cts_sum_min_thresh"]
        cts_sum_max_thresh=parameters["cts_sum_max_thresh"]
        max_jitter=args.max_jitter
        negative_sampling_ratio=args.negative_sampling_ratio
        add_revcomp=True
    elif mode=="valid":
        # no jitter at valid time - we are testing only at summits
        # no reverse complementation at valid time
        # fixed negative test for validation
        nonpeak_regions=nonpeak_regions.sample(frac=args.negative_sampling_ratio, replace=False, random_state=args.seed)
        cts_sum_min_thresh=parameters["cts_sum_min_thresh"]
        cts_sum_max_thresh=parameters["cts_sum_max_thresh"]
        max_jitter=0
        add_revcomp=False
        negative_sampling_ratio=1.0 # already subsampled
    elif mode=="test":
        # no jitter at valid time - we are testing only at summits
        # no reverse complementation at valid time
        # no subsampling of negatives - test on all positives and negatives
        # no filtering of data points
        max_jitter=0
        add_revcomp=False
        negative_sampling_ratio=1.0
        cts_sum_min_thresh="None" 
        cts_sum_max_thresh="None"
    else:
        print("mode not defined - only train, test, valid are allowed")
    return cts_sum_min_thresh,cts_sum_max_thresh,max_jitter,add_revcomp,nonpeak_regions,negative_sampling_ratio

def initialize_generators(args, mode, generator_name, parameters, return_coords):

    assert(args.inputlen%2==0)
    assert(args.outputlen%2==0)

    #defaults
    peak_regions=None
    nonpeak_regions=None

    # get only those peak/non peak regions corresponding to train/valid/test set
    if args.peaks != "None":
        print("loading peaks...")
        peak_regions=pd.read_csv(args.peaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
        peak_regions, chroms=splits.get_bed_regions_for_fold_split(peak_regions, args.fold, mode)

    if args.nonpeaks != "None":
        print("loading nonpeaks...")
        nonpeak_regions=pd.read_csv(args.nonpeaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
        nonpeak_regions, chroms=splits.get_bed_regions_for_fold_split(nonpeak_regions, args.fold, mode) 

    cts_sum_min_thresh,cts_sum_max_thresh,max_jitter,add_revcomp,nonpeak_regions, negative_sampling_ratio =  fetch_values_based_on_mode(mode, args, parameters, nonpeak_regions)
  
    if generator_name=="tiledb":
        generator=tiledb_generator.TiledbGenerator(
                                        peak_regions=peak_regions,
                                        nonpeak_regions=nonpeak_regions,
                                        genome_fasta=args.genome,
                                        batch_size=args.batch_size,
                                        input_len=args.inputlen,      
                                        output_len=args.outputlen,
                                        max_jitter=max_jitter,
                                        negative_sampling_ratio=negative_sampling_ratio,
                                        cts_sum_min_thresh=cts_sum_min_thresh,
                                        cts_sum_max_thresh=cts_sum_max_thresh,
                                        tdb_array=args.tdb_array,
                                        tdb_dataset=args.tdb_dataset,
                                        seed=args.seed,
                                        add_revcomp=add_revcomp,
                                        return_coords=return_coords
                                        )

    if generator_name=="batchgen":
        generator=batchgen_generator.ChromBPNetBatchGenerator(
                                        peak_regions=peak_regions,
                                        nonpeak_regions=nonpeak_regions,
                                        genome_fasta=args.genome,
                                        batch_size=args.batch_size,
                                        inputlen=args.inputlen,                                        
                                        outputlen=args.outputlen,
                                        max_jitter=max_jitter,
                                        negative_sampling_ratio=negative_sampling_ratio,
                                        cts_sum_min_thresh=cts_sum_min_thresh,
                                        cts_sum_max_thresh=cts_sum_max_thresh,
                                        cts_bw_file=args.bigwig,
                                        seed=args.seed,
                                        add_revcomp=add_revcomp,
                                        return_coords=return_coords
                                        )
    
    return generator
