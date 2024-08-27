import tensorflow as tf
from tensorflow import keras
import argparse
import pyBigWig
import numpy as np
import pandas as pd
import bigwig_helper
import pyfaidx
import sys
from context import data_utils as data_utils
from context import load_model_wrapper as load_model_wrapper
from context import one_hot as one_hot
import h5py

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


# need full paths!
parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS! For now, this script needs to be run from within current directory. ")
parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model h5")
parser.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=True, help="Path to chrombpnet model h5")
parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the input regions are assumed to be for [start+summit-(inputlen/2):start+summit+(inputlen/2)]. Should not be piped since it is read twice!")
parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
parser.add_argument("-o", "--out-prefix", type=str, required=True, help="Output bigwig file")
parser.add_argument("-b", "--batch-size", type=int, default=64)
parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
args = parser.parse_args()

print(args)

def get_seq(peaks_df, genome, width, regions):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []
    included = []
    pass_regions = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        if len(sequence) == width:
            vals.append(sequence)
            included.append(True)
            pass_regions.append(regions[i])
        else:
            included.append(False)
    return one_hot.dna_to_one_hot(vals), np.array(included), pass_regions

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

# load data
regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
print(regions_df.head())
gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
outputlen=1000

regions = bigwig_helper.get_regions(args.regions, outputlen) # output regions

print(len(regions))

#with pyfaidx.Fasta(args.genome) as g:
#    seqs, included, regions  = get_seq(regions_df, g, inputlen, regions)
#    print(regions_df.shape)
#    regions_df = regions_df[included].reset_index(drop=True)
#    #regions = np.array(regions)[included].tolist()
#    print(regions_df.shape)
#    print(len(regions))


chrom_h5_file = h5py.File(args.chrombpnet_model, "r")
chrom_nb_h5_file = h5py.File(args.chrombpnet_model_nb, "r")



pred_logits = chrom_h5_file["predictions"]["logits"]
pred_logcts = chrom_h5_file["predictions"]["logcounts"]

pred_logits_wo_bias = chrom_nb_h5_file["predictions"]["logits"]
pred_logcts_wo_bias = chrom_nb_h5_file["predictions"]["logcounts"]

print(pred_logits.shape)

assert(len(regions)==pred_logcts_wo_bias.shape[0])
assert(len(regions)==pred_logcts.shape[0])

pred_logits = np.squeeze(pred_logits)
pred_logits_wo_bias = np.squeeze(pred_logits_wo_bias)

bigwig_helper.write_bigwig(softmax(pred_logits_wo_bias) * (np.expand_dims(np.exp(pred_logcts_wo_bias)[:,0],axis=1)), 
                           regions, 
                           gs, 
                           args.out_prefix + "_wo_bias.bw", 
                           args.out_prefix + "_wo_bias.stats.txt", 
                           debug_chr=args.debug_chr, 
                           use_tqdm=args.tqdm)

bigwig_helper.write_bigwig(softmax(pred_logits) * (np.expand_dims(np.exp(pred_logcts)[:,0],axis=1)),
                           regions,
                           gs,
                           args.out_prefix + "_w_bias.bw",
                           args.out_prefix + "_w_bias.stats.txt",
                           debug_chr=args.debug_chr,
                           use_tqdm=args.tqdm)
