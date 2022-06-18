import os

data = open("dilation_layer_tuning_same_pad.yml").read()


template0 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
	"$bias_name": "bias.h5",
	"$dil": "13",
	"$clsmall":"hepg2",
	"$cellline": "HEPG2",
	"$fold":"0",
	"$inputlen":"2114"}


#template1 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
#	"$bias_name": "bias.h5",
#	"$dil": "8",
#	"$clsmall":"hepg2",
#	"$cellline": "HEPG2",
#	"$fold":"0",
#	"$inputlen":"2114"}

template2 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
	"$bias_name": "bias.h5",
	"$dil": "10",
	"$clsmall":"k562",
	"$cellline": "K562",
	"$bias_cl": "HEPG2",
	"$fold":"0",
	"$inputlen":"2114"}

template2 = {"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0" ,
	"$bias_name": "bias.h5",
	"$dil": "8",
	"$clsmall":"gm12878",
	"$cellline": "GM12878",
	"$bias_cl": "HEPG2",
	"$fold":"0",
	"$inputlen":"2114"}

template=template2
#change inputlen and dil and bias_name

#hyper_params = [{"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0", "$fold":"0"}]


hyper_params = [
	{"$bias_setting": "HEPG2_05.24.2022_bias_128_4_1234_0.8_fold_2", "$fold":"2"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_3", "$fold":"3"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_4", "$fold":"4"}]

hyper_params = [ 
	{"$bias_setting": "HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0", "$fold":"0"}, 
	{"$bias_setting": "HEPG2_05.20.2022_bias_128_4_1234_0.8_fold_1", "$fold":"1"},
	{"$bias_setting": "HEPG2_05.24.2022_bias_128_4_1234_0.8_fold_2", "$fold":"2"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_3", "$fold":"3"},
	{"$bias_setting": "HEPG2_05.22.2022_bias_128_4_1234_0.8_fold_4", "$fold":"4"}]


for param in hyper_params:
	for key in param:
		template[key] = param[key]

	data_n = "" + data
	for key in template:
		data_n = data_n.replace(key, template[key])


	file_name = "run_"+template["$cellline"]+"_dil_"+template["$dil"]+"_fold_"+template["$fold"]+"same_inputlen.yml"
	f = open("scripts/"+file_name, "w")
	f.write(data_n)
	f.close()

	command = "kubectl create -f " + "scripts/"+file_name
	print(command)
	os.system(command)



