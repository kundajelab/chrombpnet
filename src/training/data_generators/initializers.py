import data_generators.batchgen_generator as batchgen_generator
from utils import data_utils
import pandas as pd
import splits

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def fetch_params_based_on_mode(mode, args, parameters, nonpeak_regions):

    if mode=="train":
        # read params based on used defined values
        cts_sum_min_thresh=parameters["counts_sum_min_thresh"]
        cts_sum_max_thresh=parameters["counts_sum_max_thresh"]
        max_jitter=args.max_jitter
        negative_sampling_ratio=args.negative_sampling_ratio
        add_revcomp=True
        shuffle_at_epoch_start=True

    elif mode=="valid":
        # do not jitter at valid time - we are testing only at summits
        max_jitter=0
        # no reverse complementation at valid time
        add_revcomp=False
        # fix negative test for validation
        nonpeak_regions=nonpeak_regions.sample(frac=args.negative_sampling_ratio, replace=False, random_state=args.seed)
        negative_sampling_ratio=1.0 # already subsampled
        # filter outliers at validation - we dont want the model checkpoint to be sensitive to outliers
        cts_sum_min_thresh=parameters["counts_sum_min_thresh"]
        cts_sum_max_thresh=parameters["counts_sum_max_thresh"]
        # no need to shuffle
        shuffle_at_epoch_start=False

    elif mode=="test":
        # no jitter at valid time - we are testing only at summits
        max_jitter=0
        # no reverse complementation at valid time
        add_revcomp=False
        # no subsampling of negatives - test on all positives and negatives
        negative_sampling_ratio=1.0
        # filtering of data points based on outliers - need to discuss how to do this
        cts_sum_min_thresh="None"
        cts_sum_max_thresh="None"
        # no need to shuffle
        shuffle_at_epoch_start=False
    else:
        print("mode not defined - only train, valid, test are allowed")

    return cts_sum_min_thresh,cts_sum_max_thresh,max_jitter,add_revcomp,nonpeak_regions,negative_sampling_ratio, shuffle_at_epoch_start

def initialize_generators(args, mode, parameters, return_coords):

    assert(args.inputlen%2==0)
    assert(args.outputlen%2==0)

    # defaults
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

    cts_sum_min_thresh, cts_sum_max_thresh, max_jitter, add_revcomp, nonpeak_regions, negative_sampling_ratio, shuffle_at_epoch_start =  \
                                                                            fetch_params_based_on_mode(mode, args, parameters, nonpeak_regions)

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
                                    return_coords=return_coords,
                                    shuffle_at_epoch_start=shuffle_at_epoch_start
                                    )
    
    return generator
