import pyfaidx
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
import chrombpnet.training.utils.one_hot as one_hot

def parse_args():
    parser=argparse.ArgumentParser(description="build pwm matrix from bigwig")
    parser.add_argument("-i","--bigwig", required=True,  help="generated bigiwig file")
    parser.add_argument("-g", "--genome", required=True, help="reference genome fasta")
    parser.add_argument("-op", "--output_prefix", required=True,  help="output prefix for png file storing pwm (code will append .png suffix)")
    parser.add_argument("-cr","--chr",type=str, required=True, help="chromosome to build pwm, the name should be present in the chrom sizes file and bigwig you will provide")
    parser.add_argument("-c","--chrom_sizes",type=str, required=True, help="TSV file with chromosome name in first column and size in the second column")
    parser.add_argument("-pw","--pwm_width",type=int, default=24, required=False, help="width of pwm matrix")
    return parser.parse_args()

def get_pwm_bg(seqs, cnts, pwm_width=24):
    '''
    Arguments::
        seqs: An input array of shape L x alphabet
        cnts: An input array of shape L
    
    Returns:
        motif: PPM (Position Probability Matrix) dimensions of length x alphabet.
               Entries along the alphabet axis sum to 1.
        bg: The background base frequencies
    '''
    new_seqs = []
    for i in range(pwm_width//2,cnts.shape[0]-pwm_width//2):
        if cnts[i] > 0:
            new_seqs.append(seqs[i-pwm_width//2:i+pwm_width//2,:]*cnts[i])
    motif = np.sum(new_seqs, axis=0)
    motif = motif/np.sum(motif, axis=-1, keepdims=True)
    bg = np.sum(np.sum(new_seqs, axis=0), axis=0)
    bg = bg/sum(bg)
    return motif, bg

def main(args): 

    assert(args.pwm_width % 2 ==0)

    # access files
    hg38 = pyfaidx.Fasta(args.genome)
    bw = pyBigWig.open(args.bigwig) 

    ## find given chromosome size

    chrom_sizes_dict = {line.strip().split("\t")[0]:int(line.strip().split("\t")[1]) for line in open(args.chrom_sizes).readlines()}
    chr_size = chrom_sizes_dict[args.chr]

    # fetch values in the given chromsome and for the given chromsome region
    seq = str(hg38[args.chr][0:chr_size])
    one_hot_seq = one_hot.dna_to_one_hot(seq).squeeze()
    bigwig_vals = np.nan_to_num(bw.values(args.chr,0,chr_size ))
    print("non zero bigwig entries in the given chromosome: ", np.sum(bigwig_vals>0))

    # build pwm matrix - get PPM and background
    motif, bg = get_pwm_bg(one_hot_seq, bigwig_vals, args.pwm_width)

    # use modisco utils to plot
    figsize=(20,2)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111) 
    viz_sequence.plot_weights_given_ax(ax=ax, array=viz_sequence.ic_scale(motif, background=bg))
    plt.savefig(args.output_prefix+'.png')
    
if __name__=="__main__":
    args = parse_args()
    main(args)

