import pyfaidx
import argparse

def parse_args():
    parser=argparse.ArgumentParser(description="Bin the entire genome into bins and get gc content")
    parser.add_argument("-g","--genome", required=True, help="reference genome file")
    parser.add_argument("-o","--output-prefix", required=True, help="output BED file prefix to store the gc content of binned genome. suffix .bed will be appended by the code. If the prefix contains a directory path make sure it exists.")
    parser.add_argument("-f","--inputlen", type=int,default=2114, help="inputlen to use to make bins and find gc content")
    parser.add_argument("-s","--stride", type=int,default=1000, help="stride to use for shifting the bins")
    return parser.parse_args()

def get_genomewide_gc(genome_fa, outf, width, stride):
    """
    Get GC fraction in bins of width "width" strided by "stride".

    Main speedups come from:
    - loading chromosome string using pyfaidx
    - using the str.count function for substrings
    - avoiding redundant counting

    Redundant counting is avoided by counting in bins of size "stride"
    at a time and caching the most recent values in cache. As an example:

    For width 2114 and stride 1000, when considering the 3000-4000 bin, 
    with already cached counts in 1000-2000 and 2000-3000, count in
    3000-3114 and write 1000-3114. Then count in 3114-4000, now cache 
    3000-4000 and delete 1000-2000. And repeat.
    """

    f = pyfaidx.Fasta(genome_fa, as_raw=True, rebuild=False)
    outf = open(outf, 'w')

    div = width//stride
    rem = width%stride
    stride_x_div = div * stride

    # cache will store the GC counts in the most recent
    # div bins of length stride each
    cache = [0]*div

    for chrm in f.keys():
        s = f[chrm][:].upper()

        runsum = 0
        # fill first div values
        for i in range(0, stride_x_div, stride):
            c = s.count("C", i, i+stride) + s.count("G", i, i+stride)
            cache[i//stride % div] = c
            runsum += c
        
        for i in range(div*stride, len(s)-rem, stride):
            # invariant: runsum = sum(cache)
            left_ct=0
            if rem!=0:
                left_ct = s.count("C", i, i+rem) + s.count("G", i, i+rem)
                runsum += left_ct

            outf.write("{}\t{}\t{}\t{}\n".format(chrm, i - stride_x_div, i - stride_x_div + width, round(runsum/width,2)))

            if div == 0: # stride > width, do no more
                runsum = 0
            else:
                runsum -= cache[i//stride % div]

                right_ct = s.count("C", i+rem, i+stride) + s.count("G", i+rem, i+stride)
                runsum += right_ct
                cache[i//stride % div] = left_ct+right_ct 

    f.close()
    outf.close()
    
def main():
    args = parse_args()
    get_genomewide_gc(args.genome, args.output_prefix+".bed", args.inputlen, args.stride)

    
if __name__=="__main__":
    main()
