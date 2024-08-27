import pandas as pd

tsvf = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/script/"

annotation = pd.read_csv("atac_to_chip_1.txt",sep="\t", header=None)
anng = annotation.groupby(by=[0])

for cellty in ["GM12878", "K562", "HEPG2", "H1ESC"]:
	annt = anng.get_group(cellty)		
	print(annt.head())
	ttf = tsvf+cellty.lower()+"_chipseq.tsv"
	ppf = pd.read_csv(ttf, sep=',', header=None)

	for i,r in annt.iterrows():
		vals = list(set(ppf[ppf[1].isin(r[1].split("/"))][0]))
		if len(vals)==0 and r[1]=="None":
			#print(r)
			continue
		else:
			#print(vals)
			",".join(vals)
