import pandas as pd
import pybedtools

encid="K562"

sel_bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/selected.regions.valid.merged.bed.gz"
seldata=pd.read_csv(sel_bed, sep='\t', header=None)
print(seldata.shape)

pred_bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/preds_upload/average_preds_with_ccre_vf/filtered.regions.bed.gz"
data=pd.read_csv(pred_bed, sep='\t', header=None)
print(data.head())
x = pybedtools.BedTool.from_dataframe(data[[0,1,2]])
t = x.sort().merge()
y=pybedtools.BedTool.from_dataframe(seldata)
x = t.intersect(y, v=True, f=1.0)
print("regions in pred not in sel")

try:
	output_bed = x.to_dataframe()
	print(output_bed.shape)
	print(output_bed.head())
	print(set(output_bed["chrom"]))
except:
	print("none**************")
	pass

print("regions in sel not in pred")

x = y.intersect(t, v=True, f=1.0)
#print(y)
#print(t)
try:
	output_bed = x.to_dataframe()
	print(output_bed.shape)
	print(output_bed.head())
	print(set(output_bed["chrom"]))
except:
	print("none**************")
	pass

pred_bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/mean_folds.inputs.bed.gz"
data=pd.read_csv(pred_bed, sep='\t', header=None)
print(data.head())
data[1] = data[1]+data[9]-500
data[2] = data[1]+1000
print(data.head())
x = pybedtools.BedTool.from_dataframe(data[[0,1,2]])
t = x.sort().merge()
y=pybedtools.BedTool.from_dataframe(seldata)
x = t.intersect(y, v=True, f=1.0)

print("regions in interpret not in sel")

try:
	output_bed = x.to_dataframe()
	print(output_bed.shape)
	print(output_bed.head())
	print(set(output_bed["chrom"]))
except:
	print("none**************")
	pass


print("regions in sel not in interpret")

x = y.intersect(t, v=True, f=1.0)
try:
	output_bed = x.to_dataframe()
	print(output_bed.shape)
	print(output_bed.head())
	print(set(output_bed["chrom"]))
except:
	print("none**************")
	pass


x = y.coverage(t)
output_bed = x.to_dataframe()
print(output_bed.head())

print(sum(output_bed["thickStart"]==1.0))
print(sum(output_bed["thickStart"]<1.0))
print(seldata.shape)
print(output_bed.shape)

