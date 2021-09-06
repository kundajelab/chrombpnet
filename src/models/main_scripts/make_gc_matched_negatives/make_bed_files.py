import pandas as pd
import csv
import argparse


parser=argparse.ArgumentParser()
parser.add_argument("--input_name")
args = parser.parse_args()

input=args.input_name
negatives_bed = []

for file_type in ["train","test"]:
    negatives = open(input+".%s.0.negatives"%file_type).read().split("\n")
    for i in range(0,len(negatives)-1,2):
        headers = negatives[i]
        headers = headers[1:]
        chromo, start,end, GC = headers.split("_")
        start = (int(start)+int(end))//2
        end = start + 2
        bed = [chromo,int(start),int(end)]
        negatives_bed += [bed]


negatives_bed = [[chromo,int(start),int(end),".",int(end)-int(start),".",".",".",".",1] for chromo, start, end in negatives_bed]
negatives_bed = pd.DataFrame(negatives_bed)
outputFile = input+".all.negatives.bed"
negatives_bed.to_csv(outputFile, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
