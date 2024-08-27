import argparse
import pandas as pd
import pyBigWig
import pybedtools
import numpy as np

parser=argparse.ArgumentParser(description="")
parser.add_argument("-ac","--atac_counts", required=True,  help="")
parser.add_argument("-ap","--atac_profile", required=True,  help="")
parser.add_argument("-dc","--dnase_counts", required=True,  help="")
parser.add_argument("-dp","--dnase_profile", required=True,  help="")
parser.add_argument("-cb","--chip_bed", required=True,  help="")
parser.add_argument("-cc","--chip_contrib", required=True,  help="")
parser.add_argument("-tb","--tobias_bed", required=True,  help="")
parser.add_argument("-o","--outf", required=True,  help="")
parser.add_argument("-bwo","--observed", required=True,  help="")
parser.add_argument("-dbwo","--dnobserved", required=True,  help="")
args = parser.parse_args()



if __name__=="__main__":

	a = pybedtools.example_bedtool(args.tobias_bed)
	#b = pybedtools.example_bedtool(args.chip_bed)
	bdf = pd.read_csv(args.chip_bed, sep='\t', header=None)
	bdf[1] = bdf[1] + bdf[9] - 500
	bdf[2] = bdf[1] + 1000
	bdf[bdf[1] < 0][1] = 0

	b = pybedtools.BedTool.from_dataframe(bdf)
	a_and_b = a.intersect(b, c=True, f=1.0)
	result = a_and_b.to_dataframe()
	print(result.head())
	result['label'] = result["blockSizes"] > 0

	tobia_bed = pd.read_csv(args.tobias_bed, sep="\t", header=None)
	bw_counts = pyBigWig.open(args.atac_counts)
	bw_profile = pyBigWig.open(args.atac_profile)
	dn_bw_counts = pyBigWig.open(args.dnase_counts)
	dn_bw_profile = pyBigWig.open(args.dnase_profile)
	bw_chip = pyBigWig.open(args.chip_contrib)
	bw_obs = pyBigWig.open(args.observed)
	dn_bw_obs = pyBigWig.open(args.dnobserved)


	chip_seq_vals = []
	atac_counts = []
	atac_profiles = []
	obs_vals = []
	dnase_counts = []
	dnase_profiles = []
	dnase_obs_vals = []

	tobia_bed = tobia_bed[result['label']].reset_index(drop=True)

	for i,r in tobia_bed.iterrows():
		val1 = np.sum(np.nan_to_num(bw_counts.values(r[0],r[1],r[2])))
		atac_counts.append(val1)

		val1 = np.sum(np.nan_to_num(bw_profile.values(r[0],r[1],r[2])))
		atac_profiles.append(val1)

		val1 = np.sum(np.nan_to_num(dn_bw_counts.values(r[0],r[1],r[2])))
		dnase_counts.append(val1)

		val1 = np.sum(np.nan_to_num(dn_bw_profile.values(r[0],r[1],r[2])))
		dnase_profiles.append(val1)

		val1 = np.sum(np.nan_to_num(bw_chip.values(r[0],r[1],r[2])))
		chip_seq_vals.append(val1)
		
		mid=int((r[1]+r[2])/2)
		val1 = np.sum(np.nan_to_num(bw_obs.values(r[0],mid-150,mid+150)))
		obs_vals.append(val1)


		mid=int((r[1]+r[2])/2)
		val1 = np.sum(np.nan_to_num(dn_bw_obs.values(r[0],mid-150,mid+150)))
		dnase_obs_vals.append(val1)

	tobia_bed["counts_contrib"] = atac_counts
	tobia_bed["profiles_contrib"] = atac_profiles
	tobia_bed["dnase_counts_contrib"] = dnase_counts
	tobia_bed["dnase_profiles_contrib"] = dnase_profiles
	tobia_bed["chip"] = chip_seq_vals
	tobia_bed["observed"] = obs_vals
	tobia_bed["dnase_observed"] = dnase_obs_vals
	#tobia_bed["label"] = result['label']


	tobia_bed.to_csv(args.outf, sep='\t', header=False, index=False)
