import pysam
import argparse

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content after binning the entire genome into bins")
    parser.add_argument("-g","--genome", help="reference genome file")
    parser.add_argument("-o","--output_bed", help="output BED file to store the gc content of binned genome")
    parser.add_argument("-f","--inputlen", type=int,default=2114, help="inputlen to use to find gc content")
    parser.add_argument("-s","--stride", type=int,default=50, help="stride to use for shifting the bins")
    return parser.parse_args()

def main(genome_path, out_path, inputlen, stride):
    with pysam.FastxFile(genome_path) as fh, open(out_path, "w") as fo:
        for entry in fh:
            chrom = entry.name
            seq = entry.sequence
            buffer = [0 for _ in range(inputlen)] # Sliding window of `inputlen` bases
            ngc = 0
            for seqind, base in enumerate(seq):
                # Maintain running count of G/C while iterating through chromosome
                # `base` is the final nucleotide in bin
                ind = seqind % inputlen
                ngc -= buffer[ind] # Remove base falling out of window
                if base == "G" or base == "C": # Add new base in window
                    ngc += 1
                    buffer[ind] = 1 
                else:
                    buffer[ind] = 0

                end = seqind + 1
                start = end - inputlen
                if (start % stride == 0) and (start >= 0):
                    frac = round(ngc / inputlen, 2)
                    fo.write(f"{chrom}\t{start}\t{end}\t{frac}\n")
  
if __name__=="__main__":
    args = parse_args()
    main(args.genome, args.output_bed, args.inputlen, args.stride)