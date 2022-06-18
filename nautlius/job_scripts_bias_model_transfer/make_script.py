import os

data = open("transfer_bias_atac.yml").read()


template0 = {
	"$dil": "8",
	"$clsmall":"hepg2",
	"$cellline": "HEPG2",
	"$fold":"0"
	}

template=template0

hyper_params = [ 
	{ "$cellline":"K562", "$clsmall":"k562"}, 
	{ "$cellline":"GM12878","$clsmall":"gm12878"},
	{ "$cellline":"IMR90","$clsmall":"imr90"},
	{ "$cellline":"H1ESC","$clsmall":"h1esc"}]


for param in hyper_params:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "run_"+template["$cellline"]+"_bias_dil_"+template["$dil"]+"_fold_"+template["$fold"]+".yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()

	command = "kubectl create -f " + "scripts/"+file_name
	print(command)
	os.system(command)



