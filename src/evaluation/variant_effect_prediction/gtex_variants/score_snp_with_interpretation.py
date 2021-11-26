import math
import tensorflow
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import kerasAC 
from kerasAC.generators.snp_generator import *
from kerasAC.tiledb_config import *
from scipy.special import softmax
from kerasAC.interpret.deepshap import * 
from kerasAC.interpret.profile_shap import * 
from kerasAC.vis import * 
from kerasAC.helpers.transform_bpnet_io import * 
import pandas as pd
import os
import pickle
from load_model import * 
from scipy.special import softmax,logit
from scipy.spatial.distance import jensenshannon
import argparse

parser=argparse.ArgumentParser(description="variant effect scoring on dsQTLS")
parser.add_argument("--model_file")
parser.add_argument("--output_path")
args = parser.parse_args()

output_path=args.output_path
model_file=args.model_file

flank=500

## load data

#ref_file="/mnt/data/male.hg19.fa"
ref_file="/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
inpath = os.path.join(output_path, "formatted.csv")

snps_old=pd.read_csv("/mnt/lab_data3/anusri/histone_expts/all_qtl_analysis/gtex_variants/LCL.variants.hg38.full.tsv",header=0,sep='\t')

snps=pd.read_csv("temp.intersect.atac.bed",names=snps_old.columns,sep='\t')
snps["Chr"]=snps["chr.38"]
snps["Pos0"]=snps["start.38"]
snps["POSTallele"]=snps["ref"]
snps["ALTallele"]=snps["alt"]
print(snps.shape)
# considering only those locations with one allele type
snps = snps[snps['POSTallele'].astype(str).map(len) < 2]
snps = snps[snps['ALTallele'].astype(str).map(len) < 2]
print(snps.shape)
snps["pvalue"] = 1-snps["PIP"]
snps["logratio"]=snps["beta_posterior"]
snps['rsid']=snps['Chr'].astype(str)+'_'+snps['Pos0'].astype(str)+'_'+snps['POSTallele'].astype(str)+'_'+snps['ALTallele'].astype('str')
print(snps.head()) 
print(snps.shape)

snps.to_csv(inpath,header=True,index=False,sep='\t')
print(snps.head())
print("total data size",snps.shape)


## load model

try:
   model=load_model_wrapper(model_hdf5=model_file)
except:
   model=load_model_wrapper(json_string=model_file+"model.0.arch", weights=model_file+"model.0.weights")

## define explainers

from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
from kerasAC.interpret.deepshap import *
from kerasAC.interpret.profile_shap import *
import tensorflow as tf
import keras

task_index=0
create_background_counts=create_background_atac
model_wrapper_for_counts=(model.input, model.outputs[1][:,task_index:task_index+1])
prof_explainer = create_explainer(model,ischip=False,task_index=0)
count_explainer=shap.DeepExplainer(model_wrapper_for_counts,data=create_background_counts,combine_mult_and_diffref=combine_mult_and_diffref_1d)
print("got count explainer")
print("got profile explainer")

#reference allele sequence generator 
ref_gen=SNPGenerator(bed_path=inpath,
                 chrom_col="Chr",
                 pos_col="Pos0",
                 allele_col="POSTallele",
                 rsid_col='rsid',
                 flank_size=1057,
                 ref_fasta=ref_file,
                 batch_size=200,
                 expand_dims=False)

#alternate allele sequence generator 
alt_gen=SNPGenerator(bed_path=inpath,
                 chrom_col="Chr",
                 pos_col="Pos0",
                 allele_col="ALTallele",
                 rsid_col='rsid',
                 flank_size=1057,
                 ref_fasta=ref_file,
                 batch_size=200,
                 expand_dims=False)


#get the reference allele predictions 
count_preds={} 
profile_preds={} 
profile_shap={}
count_shap={}
final_scores={}
snp_to_seq={} 

for i in range(len(ref_gen)):
    try:
        cur_x=ref_gen[i] 
    except:
        continue
    batch_rsids=cur_x[0] 
    batch_preds=model.predict(cur_x[1])
    batch_preds_profile=batch_preds[0]
    batch_preds_count=batch_preds[1] 
    count_explanations=count_explainer.shap_values(cur_x[1])[0]
    profile_explanations=prof_explainer(cur_x[1], None)

    for batch_index in range(len(batch_rsids)): 
        cur_rsid=batch_rsids[batch_index]
        snp_to_seq[cur_rsid]={} 
        snp_to_seq[cur_rsid]['ref']=cur_x[1][batch_index,:,:]
        cur_pred_profile=batch_preds_profile[batch_index,:,:]
        cur_pred_count=batch_preds_count[batch_index,:]
        count_preds[cur_rsid]={}
        count_preds[cur_rsid]['ref']=cur_pred_count[0]
        profile_preds[cur_rsid]={}
        profile_preds[cur_rsid]['ref']=cur_pred_profile 
        profile_shap[cur_rsid]={}
        profile_shap[cur_rsid]['ref']=np.sum(profile_explanations[batch_index,:,:]*cur_x[1][batch_index,:,:], axis=1)
        count_shap[cur_rsid]={}
        count_shap[cur_rsid]['ref']=np.sum(count_explanations[batch_index,:,:]*cur_x[1][batch_index,:,:], axis=1)


#get the alternate allele predictions 
for i in range(len(alt_gen)):
    try:
        cur_x=alt_gen[i] 
    except:
        continue
    batch_rsids=cur_x[0] 
    batch_preds=model.predict(cur_x[1]) 
    batch_preds_profile=batch_preds[0]
    batch_preds_count=batch_preds[1] 
    count_explanations=count_explainer.shap_values(cur_x[1])[0]
    profile_explanations=prof_explainer(cur_x[1], None)

    for batch_index in range(len(batch_rsids)): 
        cur_rsid=batch_rsids[batch_index]
        snp_to_seq[cur_rsid]['alt']=cur_x[1][batch_index,:,:]
        cur_pred_profile=batch_preds_profile[batch_index,:,:]
        cur_pred_count=batch_preds_count[batch_index,:]
        count_preds[cur_rsid]['alt']=cur_pred_count[0]
        profile_preds[cur_rsid]['alt']=cur_pred_profile 
        profile_shap[cur_rsid]['alt']=np.sum(profile_explanations[batch_index,:,:]*cur_x[1][batch_index,:,:], axis=1)
        count_shap[cur_rsid]['alt']=np.sum(count_explanations[batch_index,:,:]*cur_x[1][batch_index,:,:], axis=1)

        ref_preds=np.array(softmax(profile_preds[cur_rsid]['ref'][:,0],axis=0))
        alt_preds=np.array(softmax(profile_preds[cur_rsid]['alt'][:,0],axis=0))
        final_scores[cur_rsid]={}
        final_scores[cur_rsid]["Alt_Minus_Ref"]=count_preds[cur_rsid]['alt']-count_preds[cur_rsid]['ref'] # this is log diff
        final_scores[cur_rsid]["sum_logratio_pred"] =np.sum(np.abs(np.log2(ref_preds+1e-8)-np.log2(alt_preds+1e-8)))
        final_scores[cur_rsid]["JSD"] = jensenshannon(ref_preds.squeeze(),alt_preds.squeeze())


print("data size after predictions", len(profile_preds.keys()))

#convert dictionary to pandas df for easier manipulation 
#Store counts
count_preds_df=pd.DataFrame.from_dict(count_preds,orient='index')
count_preds_df["rsid"] = count_preds_df.index.values
count_preds_df.to_csv(os.path.join(output_path, "count_predictions_alt_and_ref.tsv"),header=True,index=True,sep='\t')

#Store profile preds 
pickle.dump(profile_preds, open(os.path.join(output_path,"profile_predictions.pkl"), "wb" )) 

#Store profile shap 
pickle.dump(profile_shap, open(os.path.join(output_path,"profile_shap.pkl"), "wb" )) 

#Store count shap
pickle.dump(count_shap, open(os.path.join(output_path,"count_shap.pkl"), "wb" )) 

#Store variant scores
final_scores_df=pd.DataFrame.from_dict(final_scores,orient='index')
final_scores_df["rsid"] = final_scores_df.index.values
merged=pd.merge(final_scores_df,snps,left_on=['rsid'],right_on=['rsid'])
merged.to_csv(os.path.join(output_path, "variant_scores.tsv"))
print(merged.shape)


'''
#get the JSD for profile preds 
jsd_dict={} 
logratio_dict={}
combine_dict={}

for rsid in profile_preds: 
    if flank != 500:
        ref_count = count_preds_df[count_preds_df["rsid"]==rsid]["ref"]
        alt_count = count_preds_df[count_preds_df["rsid"]==rsid]["alt"]

        if len(ref_count)==0 or len(alt_count)==0:
            continue

        ref_counts_track=np.array(softmax(profile_preds[rsid]['ref'][:,0],axis=0)*np.exp(ref_count.values[0]))
        alt_counts_track=np.array(softmax(profile_preds[rsid]['alt'][:,0],axis=0)*np.exp(alt_count.values[0]))

        #print(ref_counts_track.shape)
        #combine_dict[rsid] = np.sum(np.abs(np.log2(alt_counts_track+1e-8)-np.log2(ref_counts_track+1e-8)))
        
        # get sum of counts from the middle flank probs
        count_preds_df.loc[count_preds_df["rsid"]==rsid,"ref"] = np.log(np.sum(ref_counts_track[500-flank:500+flank])+1e-8)
        count_preds_df.loc[count_preds_df["rsid"]==rsid,"alt"] = np.log(np.sum(alt_counts_track[500-flank:500+flank])+1e-8)
    
        # get profile results
        ref_preds = ref_counts_track[500-flank:500+flank]/np.sum(ref_counts_track[500-flank:500+flank])
        alt_preds = alt_counts_track[500-flank:500+flank]/np.sum(alt_counts_track[500-flank:500+flank])   
        
        count_preds_df.loc[count_preds_df["rsid"]==rsid,'Alt_Minus_Ref']=count_preds_df.loc[count_preds_df["rsid"]==rsid,'alt']-count_preds_df.loc[count_preds_df["rsid"]==rsid,'ref'] 
    else:
        ref_preds=np.array(softmax(profile_preds[rsid]['ref'][:,0],axis=0))
        alt_preds=np.array(softmax(profile_preds[rsid]['alt'][:,0],axis=0))

    cur_sum_abs_logratio=np.sum(np.abs(np.log2(ref_preds+1e-8)-np.log2(alt_preds+1e-8)))
    logratio_dict[rsid]=cur_sum_abs_logratio
    cur_jsd=jensenshannon(ref_preds.squeeze(),alt_preds.squeeze())
    jsd_dict[rsid]=cur_jsd
    
jsd_df=pd.DataFrame.from_dict(jsd_dict,orient='index')
#combine_dict_df=pd.DataFrame.from_dict(combine_dict,orient='index')

jsd_df.head()
jsd_df=jsd_df.sort_values(by=[0])
print(jsd_df.head())
print(jsd_df.tail())

logratio_df=pd.DataFrame.from_dict(logratio_dict,orient='index')
logratio_df.head() 
#sort by value 
logratio_df=logratio_df.sort_values(by=[0])
print(logratio_df.head())
print(logratio_df.tail())

jsd_df['JSD']=jsd_df[0]
jsd_df=jsd_df.drop(columns=[0])
logratio_df['sum_logratio_pred']=logratio_df[0]
logratio_df=logratio_df.drop(columns=[0])

#merge the data frames
merged=pd.merge(jsd_df,count_preds_df,left_index=True,right_index=True)
merged=pd.merge(logratio_df,merged,left_index=True,right_index=True)
merged=pd.merge(merged,snps,left_index=True,right_on=['rsid'])
merged['-log10P']=-1*np.log10(merged['pvalue'])
print("after profile preds shape", merged.shape)

merged.to_csv(os.path.join(output_path, "all_preds.tsv"))
'''
