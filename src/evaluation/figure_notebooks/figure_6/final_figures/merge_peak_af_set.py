import os
import pandas as pd

maindir="/oak/stanford/groups/akundaje/projects/african_omics/outs/10_15_2022_per_ancestry_scores/peaks/formatted/"

beds=[]
for file in os.listdir(maindir):
	peak_bed=maindir+file
	print(peak_bed)

	bed = pd.read_csv(peak_bed, sep="\t", header=None)
	print(bed.head())
	bed = bed.sort_values(by=[8], ascending=False).head(150000)
	print(bed.head())
	bed[1]  = bed[1]  + bed[9] - 100
	bed[2]  = bed[1] + 200
	print(bed.head())
	beds.append(bed)

peaks_beds = pd.concat(beds)
print(peaks_beds.shape)

import pybedtools
x = pybedtools.BedTool.from_dataframe(peaks_beds)	
c = x.sort().merge()
y = c.to_dataframe()

print(y.shape)

y.to_csv("af_caqtls/merged.peaks.af.bed", sep="\t", header=False, index=False)

