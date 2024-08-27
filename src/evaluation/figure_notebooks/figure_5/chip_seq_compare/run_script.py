import pandas as pd
import os
import multiprocessing

data = pd.read_csv("atac_to_chip.txt", sep="\t", header=None)
print(data.head())
celline = ["GM12878", "K562", "HEPG2", "H1ESC"]

#celline = ["GM12878"]


def script1(celltype, moitf, outf, chipencid):
	command="bash script.sh "+celltype+" "+moitf+" "+outf+" "+chipencid
	print(command)
	os.system(command)

for cl in celline:
	data_sub = data[data[0]==cl].reset_index()
	count = 0

	dirs="/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/"+cl+"/"+cl+"_ATAC_TOBIAS_NEW_COUNTS/"
	folders =os.listdir(dirs)
	dictss = {}
	for filed in folders:
		if filed.startswith("neg"):
			continue
		vals = filed.split(".")
		dictss[vals[1]] = filed

	arguments_list = []
	for i,r in data_sub.iterrows():
		print(i,r[2])
		if r[2] != "None":
			count+=1
			motif = dictss["pattern_"+str(i)]
			outf = "output_jan_18/"+cl+"/"+"pattern_"+str(i)+"_"+r[1]+".bed"
			if not os.path.isfile(outf):
				script1(cl, motif, outf, r[2])
			#arguments_list.append((cl, motif, outf, r[2]))
		if count==8:
			break

		#with multiprocessing.Pool() as pool:
		#	pool.starmap(script1, arguments_list)
