import json
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_dir", type=str, required=True, help="Path to store the fold information")
args = parser.parse_args()

# mention what the outputformat should look like for further editing

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8', 'chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
splits=dict()

splits[0]={'test':['chr1'],
           'valid':['chr10','chr8']}
splits[0]['train'] = [chrom for chrom in chroms if chrom not in splits[0]['test']+splits[0]['valid']]
json.dump(splits[0], open(os.path.join(args.output_dir, "fold_0.json"),"w"), indent=4)

splits[1]={'test':['chr19','chr2'],
           'valid':['chr1']}
splits[1]['train'] = [chrom for chrom in chroms if chrom not in splits[1]['test']+splits[1]['valid']]
json.dump(splits[1], open(os.path.join(args.output_dir, "fold_1.json"),"w"), indent=4)

splits[2]={'test':['chr3','chr20'],
            'valid':['chr19','chr2']}
splits[2]['train'] = [chrom for chrom in chroms if chrom not in splits[2]['test']+splits[2]['valid']]
json.dump(splits[2], open(os.path.join(args.output_dir, "fold_2.json"),"w"), indent=4)


splits[3]={'test':['chr13','chr6','chr22'],
            'valid':['chr3','chr20']}
splits[3]['train'] = [chrom for chrom in chroms if chrom not in splits[3]['test']+splits[3]['valid']]
json.dump(splits[3], open(os.path.join(args.output_dir, "fold_3.json"),"w"), indent=4)


splits[4]={'test':['chr5','chr16','chrY'],
           'valid':['chr13','chr6','chr22']}
splits[4]['train'] = [chrom for chrom in chroms if chrom not in splits[4]['test']+splits[4]['valid']]
json.dump(splits[4], open(os.path.join(args.output_dir, "fold_4.json"),"w"), indent=4)
