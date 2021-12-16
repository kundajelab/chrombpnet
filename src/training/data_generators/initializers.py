import data_generators.batchgen_generator as batchgen_generator
from utils import data_utils
import pandas as pd
import json

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def fetch_data_and_model_params_based_on_mode(mode, args, parameters, nonpeak_regions):

    if mode=="train": 
        max_jitter=args.max_jitter
        negative_sampling_ratio=args.negative_sampling_ratio
        add_revcomp=True
        shuffle_at_epoch_start=True
        inputlen=int(parameters["inputlen"])
        outputlen=int(parameters["outputlen"])

    elif mode=="valid":
        # do not jitter at valid time - we are testing only at summits
        max_jitter=0
        # no reverse complementation at valid time
        add_revcomp=False
        # fix negative test for validation
        nonpeak_regions=nonpeak_regions.sample(frac=args.negative_sampling_ratio, replace=False, random_state=args.seed)
        negative_sampling_ratio=1.0 # already subsampled
        # no need to shuffle
        shuffle_at_epoch_start=False
        inputlen=int(parameters["inputlen"])
        outputlen=int(parameters["outputlen"])
        

    elif mode=="test":
        # no jitter at valid time - we are testing only at summits
        max_jitter=0
        # no reverse complementation at test time
        add_revcomp=False
        # no subsampling of negatives - test on all positives and negatives
        negative_sampling_ratio=1.0
        # no need to shuffle
        shuffle_at_epoch_start=False
        # read input/output length
        inputlen=args.inputlen
        outputlen=args.outputlen

    else:
        print("mode not defined - only train, valid, test are allowed")

    return inputlen, outputlen, add_revcomp, max_jitter, shuffle_at_epoch_start, nonpeak_regions, negative_sampling_ratio, 


def get_bed_regions_for_fold_split(bed_regions, mode, splits_dict):
    chroms_to_keep=splits_dict[mode]
    bed_regions_to_keep=bed_regions[bed_regions["chr"].isin(chroms_to_keep)]
    print("got split:"+str(mode)+" for bed regions:"+str(bed_regions_to_keep.shape))
    return bed_regions_to_keep, chroms_to_keep

def initialize_generators(args, mode, parameters, return_coords):

    # defaults
    peak_regions=None
    nonpeak_regions=None

    # get only those peak/non peak regions corresponding to train/valid/test set
    splits_dict=json.load(open(args.chr_fold_path))

    if args.peaks.lower() != "none":
        print("loading peaks...")
        peak_regions=pd.read_csv(args.peaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
        peak_regions, chroms=get_bed_regions_for_fold_split(peak_regions, mode, splits_dict)

    if args.nonpeaks.lower() != "none":
        print("loading nonpeaks...")
        nonpeak_regions=pd.read_csv(args.nonpeaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
        nonpeak_regions, chroms=get_bed_regions_for_fold_split(nonpeak_regions, mode, splits_dict) 

    inputlen, outputlen, \
    add_revcomp, max_jitter, shuffle_at_epoch_start, \
    nonpeak_regions, negative_sampling_ratio  =  fetch_data_and_model_params_based_on_mode(mode, args, parameters, nonpeak_regions)

    assert(inputlen%2==0)
    assert(outputlen%2==0)

    generator=batchgen_generator.ChromBPNetBatchGenerator(
                                    peak_regions=peak_regions,
                                    nonpeak_regions=nonpeak_regions,
                                    genome_fasta=args.genome,
                                    batch_size=args.batch_size,
                                    inputlen=inputlen,                                        
                                    outputlen=outputlen,
                                    max_jitter=max_jitter,
                                    negative_sampling_ratio=negative_sampling_ratio,
                                    cts_bw_file=args.bigwig,
                                    seed=args.seed,
                                    add_revcomp=add_revcomp,
                                    return_coords=return_coords,
                                    shuffle_at_epoch_start=shuffle_at_epoch_start
                                    )
    
    return generator
