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

def write_predictions_h5py(output_prefix, logits, logcts, coords):

    # open h5 file for writing predictions
    output_h5_fname = "{}_predictions.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("predictions")

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_start_dset =  [int(coords[i][1]) for i in range(num_examples)]
    coords_end_dset =  [int(coords[i][2]) for i in range(num_examples)]

    dt = h5py.special_dtype(vlen=str)

    # create the "coords" group datasets
    coords_chrom_dset = coord_group.create_dataset(
        "coords_chrom", data=np.array(coords_chrom_dset, dtype=dt),
        dtype=dt, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_start_dset", data=coords_start_dset, dtype=int, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_end_dset", data=coords_end_dset, dtype=int, compression="gzip")

    # create the "predictions" group datasets
    profs_dset = pred_group.create_dataset(
        "logits",
        data=logits,
        dtype=float, compression="gzip")
    logcounts_dset = pred_group.create_dataset(
        "logcounts", data=logcts,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()


def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)


model_chrombpnet_nb = load_model_wrapper(model_hdf5=args.chrombpnet_model_nb)
model_chrombpnet = load_model_wrapper(model_hdf5=args.chrombpnet_model)

inputlen = int(model_chrombpnet_nb.input_shape[1])
outputlen = int(model_chrombpnet_nb.output_shape[0][1])

# input and output shapes should be the same for bias model and
# chrombpnet model
assert(model_chrombpnet.input_shape[1]==inputlen)
assert(model_chrombpnet.output_shape[0][1]==outputlen)

# load data
regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
print(regions_df.head())
gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
regions = bigwig_helper.get_regions(args.regions, outputlen) # output regions

assert(regions_df.shape[0]==len(regions))

if args.debug_chr is not None:
    regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
    regions = [x for x in regions if x[0]==args.debug_chr]

with pyfaidx.Fasta(args.genome) as g:
    seqs = data_utils.get_seq(regions_df, g, inputlen)

pred_logits, pred_logcts = model_chrombpnet.predict([seqs],
                                      batch_size = args.batch_size,
                                      verbose=True)


coordinates =  [[r[0], r[1], r[2]] for r in regions]

pred_logits = np.squeeze(pred_logits)
write_predictions_h5py(args.out_prefix+"_w_bias", pred_logits, pred_logcts, coordinates)

del pred_logits
del pred_logcts

pred_logits_wo_bias, pred_logcts_wo_bias = model_chrombpnet_nb.predict([seqs],
                                      batch_size = args.batch_size,
                                      verbose=True)


pred_logits_wo_bias = np.squeeze(pred_logits_wo_bias)

write_predictions_h5py(args.out_prefix+"_wo_bias", pred_logits_wo_bias, pred_logcts_wo_bias, coordinates)

f = open(args.out_prefix+"_preds_done.txt", "w")
f.write("done")
f.close()

