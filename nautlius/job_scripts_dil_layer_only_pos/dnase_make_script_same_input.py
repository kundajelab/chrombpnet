import os

data = open("dnase_dilation_layer_tuning_same_pad.yml").read()

template2 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
	"$data_type": "DNASE_SE",
	"$bias_name": "bias.h5",
	"$dil": "10",
	"$clsmall":"gm12878",
	"$cellline": "GM12878",
	"$bias_cl": "HEPG2",
	"$fold":"0",
	"$inputlen":"2114"}

template=template2

#change inputlen and dil and bias_name

hyper_params = [ 
	{"$bias_setting": "HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0", "$fold":"0"}, 
	{"$bias_setting": "HEPG2_06.09.2022_bias_128_4_1234_0.8_fold_1", "$fold":"1"},
	{"$bias_setting": "HEPG2_06.09.2022_bias_128_4_1234_0.8_fold_2", "$fold":"2"},
	{"$bias_setting": "HEPG2_06.09.2022_bias_128_4_1234_0.8_fold_3", "$fold":"3"},
	{"$bias_setting": "HEPG2_06.09.2022_bias_128_4_1234_0.8_fold_4", "$fold":"4"}]


for param in hyper_params:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "dnase_run_"+template["$cellline"]+"_dil_"+template["$dil"]+"_fold_"+template["$fold"]+"same_inputlen.yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()

	command = "kubectl create -f " + "scripts/"+file_name
	print(command)
	os.system(command)



