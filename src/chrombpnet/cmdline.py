import os
import argparse
import subprocess

from . import SRCDIR

class PathParse(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        values = os.path.abspath(values)
        setattr(namespace, self.dest, values)

def run(*args):
    subprocess.run([str(a) for a in args], check=True, cwd=SRCDIR)

def download(args):
    run("bash", "step1_download_bams_and_peaks.sh", args.data_dir)

def bigwigs(args):
    run(
        "bash", "step2_make_bigwigs_from_bams.sh", 
        args.in_bam, args.bigwig_prefix, args.data_type, args.reference_fasta, args.chrom_sizes
    )

def bias_pwm(args):
    run(
        "python", "helpers/preprocessing/analysis/build_pwm_from_bigwig.py", 
        "-i", args.bigwig, "-g", args.genome, "-o", args.output_prefix, "-c", args.chr, "-cz", args.chrom_sizes, "-pw", args.pwm_width
    )

def splits(args):
    run("python", "helpers/make_chr_splits/splits.py", "-o", args.o)

def bins(args):
    run(
        "python", "helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py", 
        "-g", args.genome, "-c", args.chrom_sizes, "-o", args.output_bed, "-f", args.inputlen, "-s", args.stride
    )

def background(args):
    run(
        "bash", "step3_get_background_regions.sh", 
        args.fasta, args.sizes, args.blacklist, args.peaks, args.inputlen, args.gc, args.output_dir, args.fold
    )

def train_bias(args):
    run(
        "bash", "step4_train_bias_model.sh", 
        args.fasta, args.bigwig, args.peaks, args.nonpeaks, args.fold, args.bias_model, args.output_dir, args.logfile
    )

def interpret_bias(args):
    run(
        "bash", "step5_interpret_bias_model.sh", 
        args.fasta, args.regions, args.model_h5, args.output_dir
    )

def train(args):
    run(
        "bash", "step6_train_chrombpnet_model.sh", 
        args.fasta, args.bigwig, args.peaks, args.nonpeaks, args.fold, args.bias_thresh, args.output_dir, args.data_type, args.logfile
    )

def interpret(args):
    run(
        "bash", "step7_interpret_chrombpnet_model.sh", 
        args.fasta, args.regions, args.model_h5, args.output_dir
    )

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser = subparsers.add_parser("download_sample_data")
    subparser.add_argument("data_dir", action=PathParse)
    subparser.set_defaults(func=download)

    subparser = subparsers.add_parser("build_bigwigs")
    subparser.add_argument("in_bam", action=PathParse)
    subparser.add_argument("bigwig_prefix", action=PathParse)
    subparser.add_argument("data_type", choices=["ATAC_PE", "DNASE_SE", "DNASE_PE", "ATAC_SE"])
    subparser.add_argument("reference_fasta", action=PathParse)
    subparser.add_argument("chrom_sizes", action=PathParse)
    subparser.set_defaults(func=bigwigs)

    subparser = subparsers.add_parser("plot_bias_pwm")
    subparser.add_argument("-i", "--bigwig", action=PathParse, help="generated bigwig file")
    subparser.add_argument("-g", "--genome", action=PathParse, help="reference genome fasta")
    subparser.add_argument("-o", "--output", action=PathParse, help="output image for storing pwm")
    subparser.add_argument("-c","--chr", type=str, help="chromosome to build pwm, the name should be present in the chrom sizes file and bigwig you will provide")
    subparser.add_argument("-cz","--chrom_sizes",type=str, action=PathParse, help="TSV file with chromosome name in first column and size in the second column")
    subparser.add_argument("-pw","--pwm_width",type=int, default=24, help="width of pwm matrix")
    subparser.set_defaults(func=bias_pwm)

    subparser = subparsers.add_parser("get_human_splits")
    subparser.add_argument("-o", action=PathParse)
    subparser.set_defaults(func=splits)

    subparser = subparsers.add_parser("bin_genome")
    subparser.add_argument("-g","--genome", action=PathParse, help="reference genome file")
    subparser.add_argument("-c","--chrom_sizes", action=PathParse, help="chromosome sizes file for reference genome (contains chr and chrom size seperated by tab)")
    subparser.add_argument("-o","--output_bed", action=PathParse, help="output BED file to store the gc content of binned genome")
    subparser.add_argument("-f","--inputlen", type=int, default=2114, help="inputlen to use to find gc content")
    subparser.add_argument("-s","--stride", type=int, default=50, help="stride to use for shifting the bins")
    subparser.set_defaults(func=bins)

    subparser = subparsers.add_parser("get_background")
    subparser.add_argument("fasta", action=PathParse)
    subparser.add_argument("sizes", action=PathParse)
    subparser.add_argument("blacklist", action=PathParse)
    subparser.add_argument("peaks", action=PathParse)
    subparser.add_argument("inputlen", type=int)
    subparser.add_argument("gc", action=PathParse)
    subparser.add_argument("output_dir", action=PathParse)
    subparser.add_argument("fold", action=PathParse)
    subparser.set_defaults(func=background)

    subparser = subparsers.add_parser("train_bias")
    subparser.add_argument("fasta", action=PathParse)
    subparser.add_argument("bigwig", action=PathParse)
    subparser.add_argument("peaks", action=PathParse)
    subparser.add_argument("nonpeaks", action=PathParse)
    subparser.add_argument("fold", action=PathParse)
    subparser.add_argument("bias_thresh", type=float)
    subparser.add_argument("output_dir", action=PathParse)
    subparser.add_argument("logfile", action=PathParse)
    subparser.set_defaults(func=train_bias)

    subparser = subparsers.add_parser("interpret_bias")
    subparser.add_argument("fasta", action=PathParse)
    subparser.add_argument("regions", action=PathParse)
    subparser.add_argument("model_h5", action=PathParse)
    subparser.add_argument("output_dir", action=PathParse)
    subparser.set_defaults(func=interpret_bias)
    
    subparser = subparsers.add_parser("train_chrombpnet")
    subparser.add_argument("fasta", action=PathParse)
    subparser.add_argument("bigwig", action=PathParse)
    subparser.add_argument("peaks", action=PathParse)
    subparser.add_argument("nonpeaks", action=PathParse)
    subparser.add_argument("fold", action=PathParse)
    subparser.add_argument("bias_model", action=PathParse)
    subparser.add_argument("output_dir", action=PathParse)
    subparser.add_argument("data_type", choices=["ATAC_PE", "DNASE_SE", "DNASE_PE", "ATAC_SE"])
    subparser.add_argument("logfile", action=PathParse)
    subparser.set_defaults(func=train)

    subparser = subparsers.add_parser("interpret_chrombpnet")
    subparser.add_argument("fasta", action=PathParse)
    subparser.add_argument("regions", action=PathParse)
    subparser.add_argument("model_h5", action=PathParse)
    subparser.add_argument("output_dir", action=PathParse)
    subparser.set_defaults(func=interpret)

    args = parser.parse_args()
    args.func(args)




