import json
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse


def read_args():
	parser = argparse.ArgumentParser(description="Make summary reports")
	parser.add_argument('-id', '--input-dir', required=True, type=str, help="directory name output by command chrombpnet bias pipeline")
	parser.add_argument('-d', '--data-type', required=True, type=str, help="assay data type - ATAC or DNASE")
	parser.add_argument('-fp', '--file-prefix', required=False, default=None, type=str, help="File prefix for output to use. All the files will be prefixed with this string if provided.")
	parser.add_argument('-icmd', '--command', required=False, default="pipeline", type=str, choices=['pipeline', 'train', "qc"], help="Choices on what to include in the report - entire pipleline, only train, only qc")

	args = parser.parse_args()
	return args
	
def train_report(fpx,prefix):

	# preprocessing defaults
	pre_hed = 'Preprocessing report'
	pre_text = 'The image should closely represent a Tn5 or DNase enzyme motif (indicates correct shift).'

	bias_model_perf_hed = 'Bias model performance in peaks'
	bias_model_perf_text = 'The pearsonr in peaks should be greater than -0.3. Otherwise the bias model has AT rich bias.'

	data = json.load(open(os.path.join(prefix,"evaluation/bias_metrics.json")))
	df = pd.json_normalize(data['counts_metrics']).round(2)
	df = pd.json_normalize(data['counts_metrics'])
	df.index = ['counts_metrics']


	df1 = pd.json_normalize(data['profile_metrics']).round(2)
	df1 = pd.json_normalize(data['profile_metrics'])
	df1.index = ['profile_metrics']

	bias_image_loc=os.path.join("./","{}bw_shift_qc.png".format(fpx))
	print(bias_image_loc)
	## training images

	train_hed = 'Training report'

	loss = pd.read_csv(os.path.join(prefix,"logs/{}chrombpnet.log".format(fpx)), sep=",", header=0)
	
	val_loss = loss["val_loss"]
	train_loss = loss["loss"]
	epochs = loss["epoch"]
	
	plt.rcParams["figure.figsize"]=4,4
	plt.figure()
	plt.plot(epochs, val_loss, label="val loss")
	plt.plot(epochs, train_loss, label="train loss")
	plt.legend(loc='best')
	plt.xlabel('Epochs')
	plt.ylabel('Total loss')
	plt.tight_layout()
	plt.savefig(os.path.join(prefix,"evaluation/{}epoch_loss.png".format(fpx)),format='png',dpi=300)
	loss_image_loc=os.path.join("./",'{}epoch_loss.png'.format(fpx))

	# 2. Combine them together using a long f-string
	html = f'''
			<body style="font-size:20px;">
				<h3>{pre_hed}</h3>
				<p>{pre_text}</p>
				<img src={{bias_image}} style="max-width: 80%;">
				<h3>{bias_model_perf_hed}</h3>
				<p>{bias_model_perf_text}</p>
				{df.to_html(classes='mystyle')}
				{df1.to_html(classes='mystyle')}
			</body>
			<body style="font-size:20px;">
				<h3>{train_hed}</h3>
				<img src={{loss_image_loc}} class="center ;"> 
			</body>
		'''	
	return html.format(bias_image=bias_image_loc,loss_image_loc=loss_image_loc)
	

def qc_report(fpx,prefix,data_type):

	## Bias factorized ChromBPNet training performance

	chrombpnet_model_perf_hed = 'ChromBPNet model performance in peaks'
	chrombpnet_model_perf_text = 'The pearsonr in peaks should be greater than 0 (higher the better). Median JSD lower the better. Median Norm JSD higher the better. '

	data = json.load(open(os.path.join(prefix,"evaluation/{}chrombpnet_metrics.json".format(fpx))))
	pdf = pd.json_normalize(data['counts_metrics']).round(2)
	pdf = pd.json_normalize(data['counts_metrics'])
	pdf.index = ['counts_metrics']

	pdf1 = pd.json_normalize(data['profile_metrics']).round(2)
	pdf1 = pd.json_normalize(data['profile_metrics'])
	pdf1.index = ['profile_metrics']
		
	## Marginal footprinting on enzyme bias
	
	marg_hed = 'ChromBPNet marginal footprints on tn5 motifs'
	marg_text = 'The peak of the profiles should be below 0.003 (indicates corrected for bias) '
	data = open(os.path.join(prefix,"evaluation/{}chrombpnet_nobias_max_bias_resonse.txt".format(fpx))).read()
	vals = data.split("_")
	marg_text1 = "The average of the peaks is: "+vals[1]+" And the model is "+vals[0]

	## TFModisco motifs learnt from ChromBPNet after bias correction (chrombpnet_nobias.h5) model 
	tf_hed = "TFModisco motifs learnt from ChromBPNet after bias correction (chrombpnet_nobias.h5) model"

	tf_text_profile = "TFModisco on Profile head - Only TF motifs should be present and no bias motifs. cwm_fwd, cwm_rev should be free from any bias motifs. The motifs top matches in TOMTOM are shown (match_0, match_1, match_2)"
	#tf_text_counts = "TFModisco on Counts head. cwm_fwd, cwm_rev should have only TF motifs.  The motifs top matches in TOMTOM are shown (match_0, match_1, match_2)"

	table_profile = open(os.path.join(prefix,"evaluation/modisco_profile/motifs.html")).read().replace("./","./modisco_profile/").replace("width=\"240\"","class=\"cover\"").replace("border=\"1\" class=\"dataframe\"","").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs")
	#table_counts = open(os.path.join(prefix,"auxiliary/interpret_subsample/modisco_counts/motifs.html")).read().replace("./","./modisco_counts/").replace("width=\"240\"","class=\"cover\"").replace("border=\"1\" class=\"dataframe\"","").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs")

	if data_type == "ATAC":
		tn5_1 = os.path.join("./","{}chrombpnet_nobias.tn5_1.footprint.png".format(fpx))
		tn5_2 = os.path.join("./","{}chrombpnet_nobias.tn5_2.footprint.png".format(fpx))
		tn5_3 = os.path.join("./","{}chrombpnet_nobias.tn5_3.footprint.png".format(fpx))
		tn5_4 = os.path.join("./","{}chrombpnet_nobias.tn5_4.footprint.png".format(fpx))
		tn5_5 = os.path.join("./","{}chrombpnet_nobias.tn5_5.footprint.png".format(fpx))

		html_table = f'''	
					<body>	
						<h3>{marg_hed}</h3>
						<p>{marg_text}</p>
						<p>{marg_text1}</p>	
						<table>
						 <thead>
							<tr>		
								<td>tn5 motif 1</td>
								<td>tn5 motif 2</td>
								<td>tn5 motif 3</td>
								<td>tn5 motif 4</td>
								<td>tn5 motif 5</td>
							</tr>
						  </thead>
						<tbody>	
							<tr>		
								<td><img src={{tn5_1}} class="cover"></td>
								<td><img src={{tn5_2}} class="cover"></td>
								<td><img src={{tn5_3}} class="cover"></td>
								<td><img src={{tn5_4}} class="cover"></td>
								<td><img src={{tn5_5}} class="cover"></td>
							</tr>
						</tbody>
						</table>
					</body>
		'''	
	elif data_type == "DNASE":
		dnase_1 = os.path.join("./","{}chrombpnet_nobias.dnase_1.footprint.png".format(fpx))
		dnase_2 = os.path.join("./","{}chrombpnet_nobias.dnase_2.footprint.png".format(fpx))

		html_table = f'''
					<body>	
						<h3>{marg_hed}</h3>
						<p>{marg_text}</p>
						<p>{marg_text1}</p>	
						<table>
						 <thead>
							<tr>		
								<td>dnase motif 1</td>
								<td>dnase motif 2</td>
							</tr>
						  </thead>
						<tbody>	
							<tr>		
								<td><img src={{dnase_1}} class="cover"></td>
								<td><img src={{dnase_2}} class="cover"></td>
							</tr>
						</tbody>
						</table>
					</body>	
		'''
	else:
		print("Unknown data type: "+data_type)

	html_perf = f'''
			
			<body style="font-size:20px;">
				<h3>{chrombpnet_model_perf_hed}</h3>
				<p>{chrombpnet_model_perf_text}</p>
				{pdf.to_html(classes='mystyle')}
				{pdf1.to_html(classes='mystyle')}
			</body>	
		'''
		
# 	html_motifs = f'''			
# 			<body style="font-size:20px;">
# 				<h3>{tf_hed}</h3>
# 				<p>{tf_text_profile}</p>
# 			</body>
# 			<body>
# 				 {table_profile}
# 			</body>
# 			<body style="font-size:20px;">
# 				<p>{tf_text_counts}</p>
# 			</body>
# 			<body>
# 				 {table_counts}
# 			</body>
# 		'''					
	html_motifs = f'''			
			<body style="font-size:20px;">
				<h3>{tf_hed}</h3>
				<p>{tf_text_profile}</p>
			</body>
			<body>
				 {table_profile}
			</body>
		'''
	html = html_perf+html_table+html_motifs
	return html.format(tn5_1=tn5_1,tn5_2=tn5_2,tn5_3=tn5_3,tn5_4=tn5_4,tn5_5=tn5_5)

def main(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
		
	prefix = args.input_dir
	data_type = args.data_type
	
	pd.set_option('colheader_justify', 'center')   # FOR TABLE <th>

	if (args.command=="pipeline"):
		page_title_text='Bias factorized ChromBPNet training and quality check report'
		main_html = f'''
		<html>
			<head>
				<title>{page_title_text}</title>
			</head>	
			<body>
				<h3>{page_title_text}</h3>
			</body>
			'''
		main_html += train_report(fpx,prefix)
		main_html += qc_report(fpx, prefix, data_type)
		
	if (args.command=="train"):
		page_title_text='Bias factorized ChromBPNet training report'
		main_html = f'''
		<html>
			<head>
				<title>{page_title_text}</title>
			</head>	
			<body>
				<h3>{page_title_text}</h3>
			</body>
			'''
		main_html += train_report(fpx,prefix)

	if (args.command=="qc"):
		page_title_text='Bias factorized ChromBPNet quality check report'
		main_html = f'''
		<html>
			<head>
				<title>{page_title_text}</title>
			</head>	
			<body>
				<h3>{page_title_text}</h3>
			</body>
			'''
		main_html += qc_report(fpx, prefix, data_type)

	end_html = f'''
	</html>
		'''	
		
	main_html+=end_html
	
	# 3. Write the html string as an HTML file
	with open(os.path.join(prefix,"evaluation/{}overall_report.html".format(fpx)), 'w') as f:
		f.write(main_html)

	from weasyprint import HTML, CSS
	css = CSS(string='''
		@page {
                size: 450mm 700mm;
                mnargin: 0in 0in 0in 0in;
                }
        .center {
                 max-width:25%;
                 max-height:5%;
                 display: block;
                 margin-left: auto;
                 margin-right: auto;
                 }
        .cover {
    		 width: 100%;
   			 display: block;
			}
			
		table {
  			   font-size: 11pt;
 			   font-family: Arial;
 			   text-align: center;
 		       width: 100%;
 		       border-collapse: collapse;
 			   border: 1px solid silver;
			}
			
        ''')
	HTML(os.path.join(prefix,"evaluation/{}overall_report.html".format(fpx))).write_pdf(os.path.join(prefix,"evaluation/{}overall_report.pdf".format(fpx)), stylesheets=[css])

if __name__=="__main__":
	args=read_args()
	main(args)