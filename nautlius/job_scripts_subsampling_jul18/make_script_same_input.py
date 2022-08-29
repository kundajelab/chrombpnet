import os

data = open("subsampling.yml").read()


template2 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
	"$fold":"0",
	"$datatype":"ATAC_PE",
	"$biasdatatype":"ATAC_PE",
	"$biascellline":"HEPG2",
	"$bias_name": "bias.h5",
	"$seed": "1234",
	"$cellline":"GM12878_250M",
	"$clsmall": "gm12878.250m"}


template=template2



hyper_params = [ 
	{"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0", "$fold":"0"}, 
	{"$bias_setting": "HEPG2_05.20.2022_bias_128_4_1234_0.8_fold_1", "$fold":"1"},
	{"$bias_setting": "HEPG2_05.24.2022_bias_128_4_1234_0.8_fold_2", "$fold":"2"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_3", "$fold":"3"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_4", "$fold":"4"}]


for param in hyper_params[1:]:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "run_"+template["$cellline"]+"_fold_"+template["$fold"]+".yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()

	command = "kubectl create -f " + "scripts/"+file_name

	print(command)
	os.system(command)



