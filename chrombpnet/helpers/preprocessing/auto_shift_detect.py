import argparse
import pyBigWig
import pyfaidx
import subprocess
import pandas as pd
import numpy as np
import itertools
import warnings
from termcolor import colored
import os
from modisco.visualization import viz_sequence
import chrombpnet.training.utils.one_hot as one_hot
from chrombpnet.data import DefaultDataFile, get_default_data_path

def parse_args():
    parser=argparse.ArgumentParser(description="Automatically detect enzyme shift of input BAM/fragment/tagAlign File")
    parser.add_argument('-g','--genome', required=True, type=str, help="reference genome file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-ibam', '--input-bam-file', type=str, help="Input BAM file")
    group.add_argument('-ifrag', '--input-fragment-file', type=str, help="Input fragment file")
    group.add_argument('-itag', '--input-tagalign-file', type=str, help="Input tagAlign file")
    parser.add_argument('-d', '--data-type', required=True, type=str, choices=['ATAC', 'DNASE'], help="assay type")
    parser.add_argument('--ATAC-ref-path', type=str, default=None, help="Path to ATAC reference motifs (ATAC.ref.motifs.txt used by default)")
    parser.add_argument('--DNASE-ref-path', type=str, default=None, help="Path to DNASE reference motifs (DNASE.ref.motfis.txt used by default)")
    parser.add_argument('--num-samples', type=int, default=10000, help="Number of reads to sample from BAM/fragment file")
    args = parser.parse_args()
    return args

def ic_scale(pwm):
    # renormalize just in case
    pwm = pwm/np.sum(pwm, axis=-1, keepdims=True)
    return viz_sequence.ic_scale(pwm, background=[.25]*4)

def convolve(to_scan, longer_seq):
    # Convolve to_scan matrix against longer_seq matrix
    vals = []
    for i in range(len(longer_seq) - len(to_scan) + 1):
        vals.append(np.sum(to_scan*longer_seq[i:i+len(to_scan)]))
    return vals
    
def stream_filtered_tagaligns(src_tagaligns_stream, genome_file, out_stream):
    '''
    Given a tagalign subprocess stream and reference genome file, filters
    out any reads in chromosomes not included in the reference. Reads in the
    reference chromosomes are sent to the specified output stream.
    
    Returns:
        Boolean. Indicates whether any reads not in the reference fasta were 
        detected.
    '''
    has_unknown_chroms = False
    with pyfaidx.Fasta(genome_file) as g:
        for line in iter(src_tagaligns_stream.stdout.readline, b''):
            tagalign_chrom = line.decode('utf-8').strip().split('\t')[0]
            if tagalign_chrom in g.keys():
                out_stream.stdin.write(line)
            else:
                has_unknown_chroms = True
        src_tagaligns_stream.stdout.close()
    return has_unknown_chroms

def is_gz_file(filepath):
    # https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def bam_to_tagalign_stream(bam_path):
    p = subprocess.Popen(["bedtools", "bamtobed", "-i", bam_path], stdout=subprocess.PIPE)
    return p

def fragment_to_tagalign_stream(fragment_file_path):
    """
    Expected format for fragment file: tsv with columns chr, start, end and optionally more columns
    """
    frag_is_gz = is_gz_file(fragment_file_path)
    read_method = "zcat " if frag_is_gz else "cat "
    cmd = read_method + fragment_file_path + """ | awk -v OFS="\\t" '{print $1,$2,$3,1000,0,"+"; print $1,$2,$3,1000,0,"-"}'"""
    p = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    return p

def tagalign_stream(tagalign_file_path):
    """
    Expected format for tagAlign: tsv with columns chr, start, end, ignored, ignored, strand
    """
    ta_is_gz = is_gz_file(tagalign_file_path)
    p = subprocess.Popen(["zcat" if ta_is_gz else "cat", tagalign_file_path], stdout=subprocess.PIPE)
    return p

def sample_reads(bam_path, fragment_file_path, tagalign_file_path, num_samples, genome_fasta_path):
    # only one of bam, fragment, tagalign is not None
    if bam_path:
        p1 = bam_to_tagalign_stream(bam_path)
    elif fragment_file_path:
        p1 = fragment_to_tagalign_stream(fragment_file_path)
    elif tagalign_file_path:
        p1 = tagalign_stream(tagalign_file_path)

    # num_samples is per strand, so multiply by 2
    p2 = subprocess.Popen(["shuf", "-n", str(2*num_samples)], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    has_unknown_chroms = stream_filtered_tagaligns(p1, genome_fasta_path, p2)
    output = p2.communicate()[0]
    
    if has_unknown_chroms:
        msg = '!!! WARNING: Input reads contain chromosomes not in the reference' \
            ' genome fasta provided. Please ensure you are using the correct' \
            ' reference genome. If you are confident you are using the correct reference' \
            ' genome, you can safely ignore this message.'
        warnings.warn(colored(msg, 'red'))

    output = [x.split("\t") for x in output.decode('utf-8').split("\n")[:-1]]
    reads = pd.DataFrame(output)
    reads = reads.rename({0:'chr', 1:'start', 2:'end', 3:'x1', 4:'x2', 5:'strand'}, axis=1)

    plus_reads = reads[reads['strand']=="+"].iloc[:,:3]
    minus_reads = reads[reads['strand']=="-"].iloc[:,:3]

    return plus_reads, minus_reads


def get_pwms(plus_reads, minus_reads, genome_file):
    plus_seqs = []
    minus_seqs = []
    with pyfaidx.Fasta(genome_file) as g:
        for _, x in plus_reads.iterrows():
            cur = str(g[x['chr']][int(x['start'])-20:int(x['start'])+20])
            if len(cur)==40: # observed edge cases in non-canonical chr e.g. chrEBV
                plus_seqs.append(cur)
        for _, x in minus_reads.iterrows():
            cur = str(g[x['chr']][int(x['end'])-20:int(x['end'])+20])
            if len(cur)==40:
                minus_seqs.append(cur)
    
    plus_pwm = one_hot.dna_to_one_hot(plus_seqs).mean(0)
    plus_pwm = (plus_pwm/np.sum(plus_pwm, axis=-1, keepdims=True))
    minus_pwm = one_hot.dna_to_one_hot(minus_seqs).mean(0)
    minus_pwm = (minus_pwm/np.sum(minus_pwm, axis=-1, keepdims=True))

    return plus_pwm, minus_pwm

def get_ref_pwms(ref_motifs_path):
    """
    Expected format of file: 
    Contains motifs one position per line, 4 columns per base tab-separated. First 
    line of each motif starts with ">" followed by name. "_" is used to name motifs 
    and name should end with "_plus" or "_minus".

    For ATAC motifs the reference motifs were constructed using the get_pwms 
    function on +4/-5 shifted tagalign files and then taking central 20 bases.
    """
    pwms = {'+':{}, '-':{}}
    cur_orient = None
    cur_motif = None
    with open(ref_motifs_path) as f:
        for x in f:
            x = x.strip()
            if x.startswith(">"):
                # format current as numpy array before starting new
                if cur_motif is not None:
                    pwms[cur_orient][cur_motif] = np.array(pwms[cur_orient][cur_motif])

                if x.endswith("_plus"):
                    cur_orient = "+"
                elif x.endswith("_minus"):
                    cur_orient = "-"
                else:
                    raise ValueError("Invalid reference motif file")
                cur_motif = x[1:]
                pwms[cur_orient][cur_motif] = []
            else:
                pwms[cur_orient][cur_motif].append([float(y) for y in x.split('\t')])
    pwms[cur_orient][cur_motif] = np.array(pwms[cur_orient][cur_motif])

    return pwms['+'], pwms['-']

def compute_shift_ATAC(ref_plus_pwms, ref_minus_pwms, plus_pwm, minus_pwm):
    plus_shifts = set()
    minus_shifts = set()

    for x in ref_plus_pwms:
        # 14 is the value when comparing unshifted BAM pwm
        shift = 14 - np.argmax(convolve(ic_scale(ref_plus_pwms[x]), ic_scale(plus_pwm)))
        plus_shifts.add(shift)
    for x in ref_minus_pwms:
        shift = 5 - np.argmax(convolve(ic_scale(ref_minus_pwms[x]), ic_scale(minus_pwm)))
        minus_shifts.add(shift)

    if len(plus_shifts) != 1 or len(minus_shifts) != 1:
        raise ValueError("Input file shifts inconsistent. Please post an issue")
    
    plus_shift = list(plus_shifts)[0]
    minus_shift = list(minus_shifts)[0]

    if (plus_shift,minus_shift) not in [(0,0)]+ list(itertools.product([3,4,5],[-4,-5,-6])):
        raise ValueError("Input shift is non-standard ({:+}/{:+}). Please post an Issue.".format(plus_shift, minus_shift))

    return plus_shift, minus_shift

def compute_shift_DNASE(ref_plus_pwms, ref_minus_pwms, plus_pwm, minus_pwm):
    plus_shifts = set()
    minus_shifts = set()

    for x in ref_plus_pwms:
        # 10 is the value when comparing unshifted BAM pwm
        shift = 10 - np.argmax(convolve(ic_scale(ref_plus_pwms[x]), ic_scale(plus_pwm)))
        plus_shifts.add(shift)
    for x in ref_minus_pwms:
        shift = 10 - np.argmax(convolve(ic_scale(ref_minus_pwms[x]), ic_scale(minus_pwm)))
        minus_shifts.add(shift)

    if len(plus_shifts) != 1 or len(minus_shifts) != 1:
        raise ValueError("Input file shifts inconsistent. Please post an Issue")

    plus_shift = list(plus_shifts)[0]
    minus_shift = list(minus_shifts)[0]

    if (plus_shift,minus_shift) not in [(0,0), (0,1)]:
        raise ValueError("Input shift is non-standard ({:+}/{:+}). Please post an Issue.".format(plus_shift, minus_shift))

    return plus_shift, minus_shift 


def compute_shift(input_bam_file, input_fragment_file, input_tagalign_file, num_samples, genome_fasta_path, data_type, ref_motifs_file):
    # only one of the 3 inputs should be non None
    assert (input_bam_file is None) + (input_fragment_file is None) + (input_tagalign_file is None) == 2, "Only one input file!"

    sampled_plus_reads, sampled_minus_reads = sample_reads(input_bam_file, input_fragment_file, input_tagalign_file, num_samples, genome_fasta_path)

    plus_pwm, minus_pwm = get_pwms(sampled_plus_reads, sampled_minus_reads, genome_fasta_path)
    ref_plus_pwms, ref_minus_pwms = get_ref_pwms(ref_motifs_file)

    if data_type=="ATAC":
        plus_shift, minus_shift = compute_shift_ATAC(ref_plus_pwms, ref_minus_pwms, plus_pwm, minus_pwm)
    elif data_type=="DNASE":
        plus_shift, minus_shift = compute_shift_DNASE(ref_plus_pwms, ref_minus_pwms, plus_pwm, minus_pwm)

    return plus_shift, minus_shift


def main():
    args = parse_args()

    if args.data_type=="ATAC":
        ref_motifs_file = args.ATAC_ref_path
        if ref_motifs_file is None:
            ref_motifs_file = get_default_data_path(DefaultDataFile.atac_ref_motifs)
    elif args.data_type=="DNASE":
        ref_motifs_file = args.DNASE_ref_path
        if ref_motifs_file is None:
            ref_motifs_file = get_default_data_path(DefaultDataFile.dnase_ref_motifs) 
    plus_shift, minus_shift = compute_shift(args.input_bam_file, 
            args.input_fragment_file, 
            args.input_tagalign_file,
            args.num_samples,
            args.genome,
            args.data_type,
            ref_motifs_file)
    

if __name__=="__main__":
    main()
