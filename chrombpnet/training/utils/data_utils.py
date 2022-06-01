import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from chrombpnet.training.utils import one_hot


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

def get_seq_cts_coords(peaks_df, genome, bw, input_width, output_width, peaks_bool):

    seq = get_seq(peaks_df, genome, input_width)
    cts = get_cts(peaks_df, bw, output_width)
    coords = get_coords(peaks_df, peaks_bool)
    return seq, cts, coords

def load_data(bed_regions, nonpeak_regions, genome_fasta, cts_bw_file, inputlen, outputlen, max_jitter):
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
    
    if nonpeak_regions is not None:
        train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords = get_seq_cts_coords(nonpeak_regions,
                                              genome,
                                              cts_bw,
                                              inputlen,
                                              outputlen,
                                              peaks_bool=0)



    cts_bw.close()
    genome.close()

    return (train_peaks_seqs, train_peaks_cts, train_peaks_coords,
            train_nonpeaks_seqs, train_nonpeaks_cts, train_nonpeaks_coords)
