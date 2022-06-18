import os

data = open("transfer_bias_dnase.yml").read()


template0 = {
	"$dil": "8",
	"$clsmall":"hepg2",
	"$cellline": "HEPG2",
	"$fold":"0"
	}

template=template0

hyper_params = [ 
	{ "$cellline":"K562", "$clsmall":"k562", "$data_type":"DNASE_PE"}, 
	{ "$cellline":"GM12878","$clsmall":"gm12878", "$data_type":"DNASE_SE"},
	{ "$cellline":"IMR90","$clsmall":"imr90", "$data_type":"DNASE_SE"},
	{ "$cellline":"H1ESC","$clsmall":"h1esc", "$data_type":"DNASE_SE"}]


for param in hyper_params:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "run_dnase"+template["$cellline"]+"_bias_dil_"+template["$dil"]+"_fold_"+template["$fold"]+".yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()

	command = "kubectl create -f " + "scripts/"+file_name
	print(command)
	os.system(command)



