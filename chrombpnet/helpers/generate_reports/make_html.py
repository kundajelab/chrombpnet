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
	parser.add_argument('-hp', '--html_prefix', required=False, default="./", help="The html prefix to use for the html file output.")

	args = parser.parse_args()
	return args
	
def train_report(fpx,prefix):

	# preprocessing defaults
	pre_hed = 'Preprocessing report'
	pre_text = 'The image below should look closely like a Tn5 or DNase bias enzyme motif. '

	bias_model_perf_hed = 'Bias model performance in peaks'
	bias_model_perf_text = '<b> Counts Metrics: </b> \
	The pearsonr in peaks should be greater than -0.3 (otherwise the bias model could potentially be capturing AT bias).  \
	MSE (Mean Squared Error) will be high in peaks. <br> <br> <b> Profile Metrics:</b>  Median JSD (Jensen Shannon Divergence between observed and predicted) lower the better. \
	Median norm JSD is median of the min-max normalized JSD where min JSD is the worst case JSD i.e JSD of observed with uniform profile and \
	max JSD is the best case JSD i.e 0. Median norm JSD is higher the better. Both JSD and median norm JSD are sensitive to read-depth. \
	Higher read-depth results in better metrics. \
	<br> \
	<br> \
	<b> What to do if your pearsonr in peaks is less than -0.3? </b> \
	In the range of -0.3 to -0.5 please be wary of your chrombpnet_wo_bias.h5 TFModisco results showing lots of GC rich motifs (> 3 in the top-10). If this is not the case \
	you can continue using the chrombpnet_wo_bias.h5. If you end up seeing a lot of GC rich motifs it is likely that bias model has learnt a different GC distribution than your GC-content in peaks. \
	If you are transferring a bias model from a different sample you can consider using a different bias model or <a href=\"https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training\">training a bias model</a> for this sample. \
	If you have trained a bias model for this sample and encounter this you might have to  \
	increase the bias_threshold_factor argument input to the <em>chrombpnet bias pipeline</em> or <em>chrombpnet bias train</em> command used in training the bias model \
	and retrain a new bias model. For more intuition about this argument refer to the <a href=\"https://github.com/kundajelab/chrombpnet/wiki/FAQ\">FAQ</a> section in wiki. \
	If the value is less than -0.5 the pipline will automatically throw an error.'
	
	data = json.load(open(os.path.join(prefix,"evaluation/bias_metrics.json")))
	df = pd.json_normalize(data['counts_metrics']).round(2)
	df = pd.json_normalize(data['counts_metrics'])
	df.index = ['counts_metrics']
	df = df[df.columns.drop(list(df.filter(regex='peaks_and_nonpeaks')))]
	df = df[df.columns.drop(list(df.filter(regex='spearmanr')))]

	df1 = pd.json_normalize(data['profile_metrics']).round(2)
	df1 = pd.json_normalize(data['profile_metrics'])
	df1.index = ['profile_metrics']
	df1 = df1[df1.columns.drop(list(df1.filter(regex='peaks_and_nonpeaks')))]
	df1 = df1[df1.columns.drop(list(df1.filter(regex='spearmanr')))]

	bias_image_loc=os.path.join("./","{}bw_shift_qc.png".format(fpx))

	## training images

	train_hed = 'Training report'
	train_text = 'The val loss (validation loss) will decrease and saturate after a few epochs.'

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
				<img src={{bias_image}} style="max-width: 60%;">
				<h3>{bias_model_perf_hed}</h3>
				<p>{bias_model_perf_text}</p>
				{df.to_html(classes='mystyle')}
				{df1.to_html(classes='mystyle')}
			</body>
			<body style="font-size:20px;">
				<h3>{train_hed}</h3>
				<p>{train_text}</p>
				<img src={{loss_image_loc}} class="center", style="max-width: 20%">
			</body>

		'''	
	return html.format(bias_image=bias_image_loc,loss_image_loc=loss_image_loc)
	

def qc_report(fpx,prefix,data_type):

	## Bias factorized ChromBPNet training performance

	chrombpnet_model_perf_hed = 'ChromBPNet model performance in peaks'
	chrombpnet_model_perf_text = '<b> Counts Metrics: </b> The pearsonr in peaks should be greater than 0.5 (higher the better). \
	MSE (Mean Squared Error) will be low in peaks. <br> <br> <b> Profile Metrics:</b>  Median JSD (Jensen Shannon Divergence between observed and predicted) lower the better. \
	Median norm JSD is median of the min-max normalized JSD where min JSD is the worst case JSD i.e JSD of observed with uniform profile and \
	max JSD is the best case JSD i.e 0. Median norm JSD is higher the better. Both JSD and median norm JSD are sensitive to read-depth. \
	Higher read-depth results in better metrics.'


	data = json.load(open(os.path.join(prefix,"evaluation/{}chrombpnet_metrics.json".format(fpx))))
	pdf = pd.json_normalize(data['counts_metrics']).round(2)
	pdf = pd.json_normalize(data['counts_metrics'])
	pdf.index = ['counts_metrics']
	pdf = pdf[pdf.columns.drop(list(pdf.filter(regex='peaks_and_nonpeaks')))]
	pdf = pdf[pdf.columns.drop(list(pdf.filter(regex='spearmanr')))]

	pdf1 = pd.json_normalize(data['profile_metrics']).round(2)
	pdf1 = pd.json_normalize(data['profile_metrics'])
	pdf1.index = ['profile_metrics']
	pdf1 = pdf1[pdf1.columns.drop(list(pdf1.filter(regex='peaks_and_nonpeaks')))]
	pdf1 = pdf1[pdf1.columns.drop(list(pdf1.filter(regex='spearmanr')))]
		
	## Marginal footprinting on enzyme bias
	
	marg_hed = 'ChromBPNet marginal footprints on tn5 motifs'
	marg_text = 'The marginal footprints are the response of the ChromBPNet no bias model to the hetergenous bias motifs. \
	If the bias correction is complete the max of the profiles below should be below 0.003 on all the \
	bias motifs. '
	data = open(os.path.join(prefix,"evaluation/{}chrombpnet_nobias_max_bias_resonse.txt".format(fpx))).read()
	vals = data.split("_")
	marg_text1 = "For your convenience we calculate here the average of the max of the profiles: "+vals[1]+" And the model according to this is <b>"+vals[0]+"</b> \
	<br> \
	<br> \
	<b> What to do if your model looks uncorrected (i.e max of profiles is greater than 0.003)? </b> <br> \
	Look at the motifs below captured by TFModisco and you should be able to see motifs that closely look like the bias motifs showing incomplete bias correction. \
	This indicates that your bias model was not completely capturing the response of the bias. We recommend that you \
	use a different pre-trained bias model. For more intuition on choosing the correct pre-trained model or retraining your bias model refer to \
	<a href=\"https://github.com/kundajelab/chrombpnet/wiki/FAQ\">FAQ</a> section in wiki."


	## TFModisco motifs learnt from ChromBPNet after bias correction (chrombpnet_nobias.h5) model 
	tf_hed = "TFModisco motifs learnt from ChromBPNet after bias correction (chrombpnet_nobias.h5) model"

	#tf_text_profile = "TFModisco on Profile head - Only TF motifs should be present and no bias motifs. cwm_fwd, cwm_rev should be free from any bias motifs. The motifs top matches in TOMTOM are shown (match_0, match_1, match_2)"
	tf_text_profile = "<b> TFModisco motifs generated from profile contribution scores of the ChromBPNet after bias correction model. </b> \
	cwm_fwd, cwm_rev are the forward and reverse complemented consolidated motifs from contribution scores in subset of random peaks. \
	These CWM motifs should be free from any bias motifs and should contain only Transcription Factor (TF) motifs.\
	For each of these motifs, we use TOMTOM to find the top-3 closest matches (match_0, match_1, match_2) from a database consisting of both \
	MEME TF motifs and heterogenous enzyme bias motifs that we have repeatedly seen in our datasets.  \
	The qvals (qval0,qval1,qval2) should be low (< 0.0001) for most of the closest TF motif hits (i.e indicating that the closest match is the correct match) - this is also generally \
	verifiable by eye as the closest match will look closely like the CWMs (atleast part of it in case of heterodimers). All the motifs in the list should look nothing like the enzyme bias motif. \
	<br> \
	<br> \
	<b> What to do if you find an obvious bias motif in the list? </b> <br> \
    This indicates that your bias model was not completely capturing the response of the bias. We recommend that you \
	use a different pre-trained bias model. For more intuition on choosing the correct pre-trained model or retraining your bias model refer to \
	<a href=\"https://github.com/kundajelab/chrombpnet/wiki/FAQ\">FAQ</a> section in wiki. \
	<br> \
	<br> \
	<b> What to do if you find an obvious bias motif in the list? </b> <br>" 

	#tf_text_counts = "TFModisco on Counts head. cwm_fwd, cwm_rev should have only TF motifs.  The motifs top matches in TOMTOM are shown (match_0, match_1, match_2)"

	table_profile = open(os.path.join(prefix,"evaluation/modisco_profile/motifs.html")).read().replace("./","./modisco_profile/").replace("width=\"240\"","width=\"240\", class=\"cover\"").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs").replace("dataframe","new")
	#table_profile = open(os.path.join(prefix,"evaluation/modisco_profile/motifs.html")).read().replace("./","./modisco_profile/").replace("width=\"240\"","class=\"cover\"").replace("border=\"1\" class=\"dataframe\"","").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs").replace("dataframe","new")
	#table_counts = open(os.path.join(prefix,"auxiliary/interpret_subsample/modisco_counts/motifs.html")).read().replace("./","./modisco_counts/").replace("width=\"240\"","class=\"cover\"").replace("border=\"1\" class=\"dataframe\"","").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs")

	html_perf = f'''
			
			<body style="font-size:20px;">
				<h3>{chrombpnet_model_perf_hed}</h3>
				<p>{chrombpnet_model_perf_text}</p>
				{pdf.to_html(classes='mystyle')}
				{pdf1.to_html(classes='mystyle')}
			</body>	
		'''
		
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
								<td><img src={{tn5_1}} class="cover" style="width: 100%; display: block;"></td>
								<td><img src={{tn5_2}} class="cover" style="width: 100%; display: block;"></td>
								<td><img src={{tn5_3}} class="cover" style="width: 100%; display: block;"></td>
								<td><img src={{tn5_4}} class="cover" style="width: 100%; display: block;"></td>
								<td><img src={{tn5_5}} class="cover" style="width: 100%; display: block;"></td>
							</tr>
						</tbody>
						</table>
					</body>
		'''	
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
								<td><img src={{dnase_1}} class="cover" style="width: 100%; display: block;"></td>
								<td><img src={{dnase_2}} class="cover" style="width: 100%; display: block;"></td>
							</tr>
						</tbody>
						</table>
					</body>	
		'''
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
		return html.format(dnase_1=dnase_1,dnase_2=dnase_2)

	else:
		print("Unknown data type: "+data_type)


		
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
			<body>
			    <style> 
  				    table.dataframe, th, td {{font-size:11pt; border:1px solid black; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
			    <style> 
  				    table.new {{font-size:11pt; border:1px solid black; width:100% ; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
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
			<body>
			    <style> 
  				    table.dataframe, th, td {{font-size:11pt; border:1px solid black; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
			    <style> 
  				    table.new {{font-size:11pt; border:1px solid black; width:100% ; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
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
			<body>
			    <style> 
  				    table.dataframe, th, td {{font-size:11pt; border:1px solid black; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
			    <style> 
  				    table.new {{font-size:11pt; border:1px solid black; width:100% ; border-collapse:collapse; text-align:center;}}
 				   th, td {{padding: 5px;}}
			    </style>
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

	with open(os.path.join(prefix,"evaluation/{}overall_report.html".format(fpx)), 'w') as f:
		f.write(main_html.replace("./",args.html_prefix))
	
if __name__=="__main__":
	args=read_args()
	main(args)
