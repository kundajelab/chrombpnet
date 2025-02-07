import pandas as pd
import os
import multiprocessing

celline = ["GM12878", "K562", "HEPG2", "H1ESC"]

#celline = ["GM12878"]


def script1(celltype, moitf, outf, chipencid):
	command="bash script_new.sh "+celltype+" "+moitf+" "+outf+" "+chipencid
	print(command)
	os.system(command)

tsvf = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/script/"

annotation = pd.read_csv("atac_to_chip_1.txt",sep="\t", header=None)
anng = annotation.groupby(by=[0])

for cellty in celline:
	cl=cellty

	count = 0
	dirs="/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/"+cl+"/"+cl+"_ATAC_TOBIAS_NEW_COUNTS/"
	folders =os.listdir(dirs)
	dictss = {}
	for filed in folders:
		if filed.startswith("neg"):
			continue
		vals = filed.split(".")
		dictss[vals[1]] = filed

	annt = anng.get_group(cellty).reset_index()
	print(annt.head())
	ttf = tsvf+cellty.lower()+"_chipseq.tsv"
	ppf = pd.read_csv(ttf, sep=',', header=None)

	for i,r in annt.iterrows():
		print(r)
		vals = list(set(ppf[ppf[1].isin(r[1].split("."))][0]))
		if len(vals)==0 and r[1]=="None":
			continue
		else:					
			count+=1
			for encid in vals:

				arguments_list = []
				motif = dictss["pattern_"+str(i)]
				os.makedirs("new_output_oct_25/"+cl+"/"+"pattern_"+str(i)+"_"+r[1]+"/", exist_ok=True)
				outf = "new_output_oct_25/"+cl+"/"+"pattern_"+str(i)+"_"+r[1]+"/"+encid+".bed"
				if not os.path.isfile(outf):
					script1(cl, motif, outf, encid)

