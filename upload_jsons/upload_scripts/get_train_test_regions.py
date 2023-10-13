import pandas as pd
import json
import argparse
import os

parser = argparse.ArgumentParser(description="do bed formatting")
parser.add_argument("-ip", "--input_peaks")
parser.add_argument("-inp", "--input_nonpeaks")
parser.add_argument("-inpt", "--input_nonpeaks_test")
parser.add_argument("-f", "--chr_fold_path")
parser.add_argument("-o","--output_path")
args = parser.parse_args()

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

peaks = pd.read_csv(args.input_peaks, sep="\t", header=None, names=NARROWPEAK_SCHEMA)
nonpeaks = pd.read_csv(args.input_nonpeaks, sep="\t", header=None, names=NARROWPEAK_SCHEMA)
nonpeaks_test =  pd.read_csv(args.input_nonpeaks_test, sep="\t", header=None, names=NARROWPEAK_SCHEMA)

splits_dict=json.load(open(args.chr_fold_path))

peaks_train=peaks[peaks["chr"].isin(splits_dict['train'])]
peaks_valid=peaks[peaks["chr"].isin(splits_dict['valid'])]
peaks_test=peaks[peaks["chr"].isin(splits_dict['test'])]

nonpeaks_train=nonpeaks[nonpeaks["chr"].isin(splits_dict['train'])]
nonpeaks_valid=nonpeaks[nonpeaks["chr"].isin(splits_dict['valid'])]

path_peaks_train = os.path.join(args.output_path,"peaks.trainingset.bed.gz")
peaks_train.to_csv(path_peaks_train,sep="\t", index=False, header=False, compression='gzip')

path_peaks_valid = os.path.join(args.output_path,"peaks.validationset.bed.gz")
peaks_valid.to_csv(path_peaks_valid,sep="\t", index=False, header=False, compression='gzip')

path_peaks_test = os.path.join(args.output_path,"peaks.testset.bed.gz")
peaks_test.to_csv(path_peaks_test,sep="\t", index=False, header=False, compression='gzip')

path_nonpeaks_train = os.path.join(args.output_path,"nonpeaks.trainingset.bed.gz")
nonpeaks_train.to_csv(path_nonpeaks_train,sep="\t", index=False, header=False, compression='gzip')

path_nonpeaks_valid = os.path.join(args.output_path,"nonpeaks.validationset.bed.gz")
nonpeaks_valid.to_csv(path_nonpeaks_valid,sep="\t", index=False, header=False, compression='gzip')

path_nonpeaks_test = os.path.join(args.output_path,"nonpeaks.testset.bed.gz")
nonpeaks_test.to_csv(path_nonpeaks_test,sep="\t", index=False, header=False, compression='gzip')


