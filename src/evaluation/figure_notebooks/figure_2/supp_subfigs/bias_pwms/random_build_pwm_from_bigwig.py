import pyfaidx
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
import one_hot

def compute_per_position_ic(ppm, background, pseudocount):
    """Compute information content at each position of ppm.
    Arguments:
        ppm: should have dimensions of length x alphabet. Entries along the
            alphabet axis should sum to 1.
        background: the background base frequencies
        pseudocount: pseudocount to be added to the probabilities of the ppm
            to prevent overflow/underflow.
    Returns:
        total information content at each positon of the ppm.
    """
    assert len(ppm.shape)==2
    assert ppm.shape[1]==len(background),\
            "Make sure the letter axis is the second axis"
    if (not np.allclose(np.sum(ppm, axis=1), 1.0, atol=1.0e-5)):
        print("WARNING: Probabilities don't sum to 1 in all the rows; this can"
              +" be caused by zero-padding. Will renormalize. PPM:\n"
              +str(ppm)
              +"\nProbability sums:\n"
              +str(np.sum(ppm, axis=1)))
        ppm = ppm/np.sum(ppm, axis=1)[:,None]

    alphabet_len = len(background)
    ic = ((np.log((ppm+pseudocount)/(1 + pseudocount*alphabet_len))/np.log(2))
          *ppm - (np.log(background)*background/np.log(2))[None,:])
    return np.sum(ic,axis=1)


def ic_scale(pwm,background):
    per_position_ic = compute_per_position_ic(
                       ppm=pwm, background=background, pseudocount=0.001)
    return pwm*(per_position_ic[:,None])

def parse_args():
    parser=argparse.ArgumentParser(description="build pwm matrix from bigwig")
    parser.add_argument("-i","--bigwig", required=True,  help="generated bigiwig file")
    parser.add_argument("-g", "--genome", required=True, help="reference genome fasta")
    parser.add_argument("-o", "--output_prefix", required=True,  help="output dir for storing pwm")
    parser.add_argument("-c","--chr",type=str, required=True, help="chromosome to build pwm, the name should be present in the chrom sizes file and bigwig you will provide")
    parser.add_argument("-cz","--chrom_sizes",type=str, required=True, help="TSV file with chromosome name in first column and size in the second column")
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

if __name__=="__main__":

    args = parse_args()

    assert(args.pwm_width % 2 ==0)

    # access files
    hg38 = pyfaidx.Fasta(args.genome)
    bw = pyBigWig.open(args.bigwig) 

    ## find given chromosome size

    chrom_sizes_dict = {line.strip().split("\t")[0]:int(line.strip().split("\t")[1]) for line in open(args.chrom_sizes).readlines()}
    chr_size = chrom_sizes_dict[args.chr]

    # fetch values in the given chromsome and for the given chromsome region
    seq = str(hg38[args.chr][0:chr_size])
    import random
    one_hot_seq = one_hot.dna_to_one_hot(seq).squeeze()
    z = np.nan_to_num(bw.values(args.chr,0,chr_size ))
    random.shuffle(z)
    bigwig_vals = z
    print("non zero bigwig entries in the given chromosome: ", np.sum(bigwig_vals>0))

    # build pwm matrix - get PPM and background
    motif, bg = get_pwm_bg(one_hot_seq, bigwig_vals, args.pwm_width)

    scaled_motif = ic_scale(motif,bg)

    np.savetxt(args.output_prefix, scaled_motif)
    # use modisco utils to plot
    #figsize=(20,2)
    #fig = plt.figure(figsize=figsize)
    #ax = fig.add_subplot(111) 
    #viz_sequence.plot_weights_given_ax(ax=ax, array=viz_sequence.ic_scale(motif, background=bg))
    #plt.savefig(args.output_prefix)
