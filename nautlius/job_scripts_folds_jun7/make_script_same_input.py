import os

data = open("all_folds_runs.yml").read()


#template0 = {
#	"$biasth":"0.4",
#	"$clsmall":"gm12878",
#	"$cellline": "GM12878",
#	"$fold":"0",
#	"$seed":"1234",
#	"$dtt": "atac",
#       "$datatype":"ATAC_PE"
#	}

template0 = {
	"$biasth":"0.3",
	"$clsmall":"imr90",
	"$cellline": "IMR90",
	"$fold":"0",
	"$seed":"1234",
	"$dtt": "atac",
        "$datatype":"ATAC_PE"
	}

#template0 = {
#	"$biasth":"0.5",
#	"$clsmall":"h1esc",
#	"$cellline": "H1ESC",
#	"$fold":"0",
#	"$seed":"1234",
#	"$dtt": "atac",
#       "$datatype":"ATAC_PE"}

#template0 = {
#	"$biasth":"0.8",
#	"$clsmall":"h1esc",
#	"$cellline": "H1ESC",
#	"$fold":"0",
#	"$seed":"1234",
#	"$dtt": "atac",
#       "$datatype":"ATAC_PE"}

#template0 = {
#	"$biasth":"0.5",
#	"$clsmall":"k562",
#	"$cellline": "K562",
#	"$fold":"0",
#	"$seed":"2356",
#	"$dtt": "atac",
#        "$datatype":"ATAC_PE"}

template=template0

#change inputlen and dil and bias_name

#hyper_params = [ 
#	{"$fold":"1"},
#	{"$fold":"2"},
#	{"$fold":"3"},
#	{"$fold":"4"}]
hyper_params = [ 
	{"$fold":"1"},
	{"$fold":"2"}]



for param in hyper_params:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "run_"+template["$cellline"]+"_fold_"+template["$fold"]+"_datatype_"+template["$datatype"]+".yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()


	command = "kubectl create -f " + "scripts/"+file_name
	print(command)
	os.system(command)



