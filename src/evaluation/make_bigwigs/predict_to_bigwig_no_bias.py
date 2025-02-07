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
import os

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


# need full paths!
parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS! For now, this script needs to be run from within current directory. ")
parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model h5")
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


model_chrombpnet = load_model_wrapper(model_hdf5=args.chrombpnet_model)


inputlen = int(model_chrombpnet.input_shape[1])
outputlen = int(model_chrombpnet.output_shape[0][1])
# input and output shapes should be the same for bias model and
# chrombpnet model
assert(model_chrombpnet.input_shape[1]==inputlen)
assert(model_chrombpnet.output_shape[0][1]==outputlen)

# load data
regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
print(regions_df.head())
gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
regions = bigwig_helper.get_regions(args.regions, outputlen) # output regions

npreg = np.array(regions)

chrms = set(npreg[:,0])
print(chrms)


if args.debug_chr is not None:
    regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
    regions = [x for x in regions if x[0]==args.debug_chr]


for chrm in chrms:

    if os.path.isfile(args.out_prefix + "_" + chrm +"_w_bias.stats.txt"):
        print(args.out_prefix + "_" + chrm +"_w_bias.stats.txt")
        continue

    chrm_regions_df = regions_df[regions_df['chr'].isin([chrm])].reset_index()
    regt = [x for x in regions if x[0]==chrm]

    print(chrm_regions_df.shape)
    print(len(regt))

    with pyfaidx.Fasta(args.genome) as g:
        seqs, included, regt  = get_seq(chrm_regions_df, g, inputlen, regt)
        chrm_regions_df = chrm_regions_df[included].reset_index(drop=True)
        #regions = np.array(regions)[included].tolist()
        print(chrm_regions_df.shape)
        print(len(regt))

    print(chrm,  len(regt), regt[0], regt[-1])

    pred_logits, pred_logcts = model_chrombpnet.predict([seqs],
                                      batch_size = args.batch_size,
                                      verbose=True)

    pred_logits = np.squeeze(pred_logits)


    bigwig_helper.write_bigwig(softmax(pred_logits) * (np.expand_dims(np.exp(pred_logcts)[:,0],axis=1)),
                           regt,
                           gs,
                           args.out_prefix + "_" + chrm + "_w_bias.bw",
                           args.out_prefix + "_" + chrm +"_w_bias.stats.txt",
                           debug_chr=args.debug_chr,
                           use_tqdm=args.tqdm)

    del pred_logits, pred_logcts, seqs
