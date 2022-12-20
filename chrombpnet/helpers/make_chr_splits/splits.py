import json
import argparse
import os
import pandas as pd

def get_parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument("-op", "--output_prefix", type=str, required=True, help="Path prefix to store the fold information (appended with .json)")
    parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="TSV file with chromosome sizes. All chromosomes from the first column of chrom sizes file are used")
    parser.add_argument("-tcr", "--test_chroms", nargs="*", type=str, required=True, help="Chromosomes to use for test")
    parser.add_argument("-vcr", "--valid_chroms", nargs="*", type=str, required=True, help="Chromosomes to use for validation")
    return parser


def main(args): 

    chrom_sizes = pd.read_csv(args.chrom_sizes,sep="\t",header=None)


    chroms = chrom_sizes[0].values
    print(args.test_chroms)

    assert (set(args.test_chroms) <= set(chroms)), "Test chromosomes are not in chrom sizes file!"
    assert (set(args.valid_chroms) <= set(chroms)), "Valid chromosomes are not in chrom sizes file!"
    assert not (bool(set(args.valid_chroms) & set(args.test_chroms))), "Test and Valid chromosomes should not share any chromosomes! "


    splits=dict()
    splits[0]={'test':args.test_chroms,
               'valid':args.valid_chroms}
    splits[0]['train'] = [chrom for chrom in chroms if chrom not in splits[0]['test']+splits[0]['valid']]
    json.dump(splits[0], open(args.output_prefix + ".json","w"), indent=4)

if __name__=="__main__":
    parser = get_parsers()
    args = parser.parse_args()    
    main(args)
    
