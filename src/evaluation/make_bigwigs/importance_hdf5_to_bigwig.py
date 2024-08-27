import argparse
import pyBigWig
import numpy as np
import deepdish
import bigwig_helper

# need full paths!
parser = argparse.ArgumentParser(description="Convert importance scores in hdf5 format to bigwig. The output can be visualised using WashU Epigenome Browser as a dynseq track. Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS!")
parser.add_argument("-h5", "--hdf5", type=str, required=True, help="HDF5 file f such that f['projected_shap']['seq'] has (N x 4 x seqlen) shape with importance score * sequence so that at each f['projected_shap']['seq'][i, :, j] has 3 zeros and 1 non-zero value")
parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the importance scores are assumed to be for [start+summit-(seqlen/2):start+summit+(seqlen/2)]")
parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output bigwig file")
parser.add_argument("-s", "--outstats", type=str, required=True, help="Output file with stats of low and high quantiles")
parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
parser.add_argument("-d", "--debug-chr", type=str, default=None, help="Run for one chromosome only (e.g. chr12) for debugging")
args = parser.parse_args()

print(args)

d = deepdish.io.load(args.hdf5, '/projected_shap/seq')

SEQLEN = d.shape[2]
assert(SEQLEN%2==0)

gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
regions = bigwig_helper.get_regions(args.regions, SEQLEN)

assert(len(regions)==d.shape[0])
bigwig_helper.write_bigwig(d.sum(1), 
                           regions, 
                           gs, 
                           args.outfile, 
                           args.outstats, 
                           debug_chr=args.debug_chr, 
                           use_tqdm=args.tqdm)
