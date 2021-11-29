import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from utils import one_hot


def get_seq(peaks_df, genome, width):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []

    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        vals.append(sequence)

    return one_hot.dna_to_one_hot(vals)


def get_cts(peaks_df, bw, width):
    """
    Fetches values from a bigwig bw, given a df with minimally
    chr, start and summit columns. Summit is relative to start.
    Retrieves values of specified width centered at summit.

    "cts" = per base counts across a region
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append(np.nan_to_num(bw.values(r['chr'], 
                                            r['start'] + r['summit'] - width//2,
                                            r['start'] + r['summit'] + width//2)))
        
    return np.array(vals)

def get_coords(peaks_df, peaks_bool):
    """
    Fetch the co-ordinates of the regions in bed file
    returns a list of tuples with (chrom, summit)
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append([r['chr'], r['start']+r['summit'], "f", peaks_bool])

    return np.array(vals)

def filter_regions(peaks_df, bw, width, peaks_bool):
    """
    Filter regions in bed file that are on edges i.e regions that cannot be used to construct
    input length + jitter length of the sequence
    """
    input_shape = peaks_df.shape[0]

    # left edge case
    filtered = np.array((peaks_df['start'] + peaks_df['summit'] - width//2) < 0)
    peaks_df = peaks_df[~filtered]
    num_filtered = sum(filtered)

    # right edge case
    chrom_to_sizes = bw.chroms()
    filtered = []
    for i, r in peaks_df.iterrows():
        if r['start'] + r['summit'] + width//2 > chrom_to_sizes[r['chr']] :
            filtered.append(True)
        else:
            filtered.append(False)
    filtered=np.array(filtered)
    peaks_df = peaks_df[~filtered]
    num_filtered += sum(filtered)

    if peaks_bool:
        print("Number of peaks input: ",input_shape)
        print("Number of peaks filtered because the input/output is on the edge: ", num_filtered)
        print("Number of peaks being used: ",peaks_df.shape[0])
    else:
        print("Number of non peaks input: ",input_shape)
        print("Number of non peaks filtered because the input/output is on the edge: ", num_filtered)
        print("Number of non peaks being used: ",peaks_df.shape[0])

    return peaks_df

def get_seq_cts_coords(peaks_df, genome, bw, input_width, output_width, peaks_bool):

    # filter regions that are on the edges of the chromosome
    peaks_df = filter_regions(peaks_df, bw, input_width, peaks_bool)

    seq = get_seq(peaks_df, genome, input_width)
    cts = get_cts(peaks_df, bw, output_width)
    coords = get_coords(peaks_df, peaks_bool)
    return seq, cts, coords

def load_data(bed_regions, nonpeak_regions, genome_fasta, cts_bw_file, inputlen, outputlen, max_jitter, cts_sum_min_thresh, cts_sum_max_thresh):
    """
    Load sequences and corresponding base resolution counts for training, 
    validation regions in peaks and nonpeaks (2 x 2 x 2 = 8 matrices).

    For training peaks/nonpeaks, values for inputlen + 2*max_jitter and outputlen + 2*max_jitter 
    are returned centered at peak summit. This allows for jittering examples by randomly
    cropping. Data of width inputlen/outputlen is returned for validation
    data.

    If outliers is not None, removes training examples with counts > outlier%ile
    """

    cts_bw = pyBigWig.open(cts_bw_file)
    genome = pyfaidx.Fasta(genome_fasta)

    train_peaks_seqs=None
    train_peaks_cts=None
    train_peaks_coords=None
    train_nonpeaks_seqs=None
    train_nonpeaks_cts=None
    train_nonpeaks_coords=None

    if bed_regions is not None:
        train_peaks_seqs, train_peaks_cts, train_peaks_coords = get_seq_cts_coords(bed_regions,
                                              genome,
                                              cts_bw,
                                              inputlen+2*max_jitter,
                                              outputlen+2*max_jitter,
                                              peaks_bool=1)

        if cts_sum_min_thresh.lower() !=  "none":
            midpoint=(outputlen+2*max_jitter)//2
            train_peaks_subset = train_peaks_cts[:,midpoint-outputlen//2:midpoint+outputlen//2].sum(-1) > float(cts_sum_min_thresh)
            train_peaks_seqs = train_peaks_seqs[train_peaks_subset]
            train_peaks_cts = train_peaks_cts[train_peaks_subset]
            train_peaks_coords = train_peaks_coords[train_peaks_subset]

        if cts_sum_max_thresh.lower() != "none":
            midpoint=(outputlen+2*max_jitter)//2
            train_peaks_subset = train_peaks_cts[:,midpoint-outputlen//2:midpoint+outputlen//2].sum(-1) < float(cts_sum_max_thresh)
            train_peaks_seqs = train_peaks_seqs[train_peaks_subset]
            train_peaks_cts = train_peaks_cts[train_peaks_subset]
            train_peaks_coords = train_peaks_coords[train_peaks_subset]

        print("Total peaks after applying the min/max thresholds: ", train_peaks_seqs.shape[0])
    
    if nonpeak_regions is not None:
        train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords = get_seq_cts_coords(nonpeak_regions,
                                              genome,
                                              cts_bw,
                                              inputlen+2*max_jitter,
                                              outputlen+2*max_jitter,
                                              peaks_bool=0)

        if cts_sum_min_thresh.lower() !=  "none":
            midpoint=(outputlen+2*max_jitter)//2
            train_nonpeaks_subset = train_nonpeaks_cts[:,midpoint-outputlen//2:outputlen+outputlen//2].sum(-1) > float(cts_sum_min_thresh)
            train_nonpeaks_seqs = train_nonpeaks_seqs[train_nonpeaks_subset]
            train_nonpeaks_cts = train_nonpeaks_cts[train_nonpeaks_subset]
            train_nonpeaks_coords = train_nonpeaks_coords[train_nonpeaks_subset]

        if cts_sum_max_thresh.lower() != "none":
            midpoint=(outputlen+2*max_jitter)//2
            train_nonpeaks_subset = train_nonpeaks_cts[:,midpoint-outputlen//2:outputlen+outputlen//2].sum(-1) < float(cts_sum_max_thresh)
            train_nonpeaks_seqs = train_nonpeaks_seqs[train_nonpeaks_subset]
            train_nonpeaks_cts = train_nonpeaks_cts[train_nonpeaks_subset]
            train_nonpeaks_coords = train_nonpeaks_coords[train_nonpeaks_subset]
        
        print("Total non peaks after applying the min/max thresholds: ", train_nonpeaks_seqs.shape[0])

    cts_bw.close()
    genome.close()

    return (train_peaks_seqs, train_peaks_cts, train_peaks_coords,
            train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords)