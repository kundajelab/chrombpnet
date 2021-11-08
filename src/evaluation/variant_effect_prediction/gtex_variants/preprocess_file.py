import pandas as pd

snp=pd.read_csv("/mnt/lab_data3/anusri/histone_expts/all_qtl_analysis/gtex_variants/LCL.variants.hg38.full.tsv",header=0,sep='\t')

index1 = snp.index[snp['ref'].astype(str).map(len) > 1].tolist()
for id in index1[0:1]:
	print("***************")
	vals = snp.loc[id,"ref"]
	uniq_id = snp.loc[id, "unique.id.hg38"]
	print(snp[snp["unique.id.hg38"]==uniq_id])
	snp.loc[id,"ref"] = vals[0]
	row = snp.loc[id].copy()
	print(snp[snp["unique.id.hg38"]==uniq_id])
	print(vals)
	for i in range(1,len(vals)):
		print(i)
		row1 = row.copy()
		row1.ref = vals[i]
		print(row1)
		snp = snp.append(row1)
		print(snp.shape)
	print(snp[snp["unique.id.hg38"]==uniq_id])
