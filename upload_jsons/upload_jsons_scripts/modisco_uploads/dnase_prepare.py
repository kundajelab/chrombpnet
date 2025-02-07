import os
import json

encids = ["IMR90_new", "H1ESC_new", "GM12878_new", "HEPG2", "K562"]

encode_id = {"HEPG2": "ENCSR149XIL",
        "K562": "ENCSR000EOT",
        "IMR90_new": "ENCSR477RTP",
        "GM12878_new": "ENCSR000EMT",
        "H1ESC_new": "ENCSR000EMU"}

ooutdir='dnase/'





for name in encids:

	encid = encode_id[name]
	args_json = {}
	args_json["experiment"] = encode_id[name]
	args_json["sequence motifs tar"] = {}

	success=True
	readme_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/dummy/chrombpnet_test/READMEs/modisco.report.README"
	if os.path.isfile(readme_file):
		args_json["sequence motifs tar"]["file.paths"] = [(readme_file, "README.md")]
		
	args_json["sequence motifs tar"]["counts"] = {"file.paths": []}
	args_json["sequence motifs tar"]["profile"] = {"file.paths": []}

	odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"
	counts_modisco=odir+name+"/merge_folds_new_may_05_24/counts/modisco_counts.h5"
	if os.path.isfile(counts_modisco):
		args_json["sequence motifs tar"]["counts"]["file.paths"].append((counts_modisco, "tfmodisco.raw_output.counts."+encid+".hd5"))
	else:
		print(counts_modisco)
		continue
	
	profile_modisco=odir+name+"/merge_folds_new_may_05_24/profile/modisco_profile.h5"
	if os.path.isfile(profile_modisco):
		args_json["sequence motifs tar"]["profile"]["file.paths"].append((profile_modisco, "tfmodisco.raw_output.profile."+encid+".hd5"))
	else:
		print(profile_modisco)
		continue
	
	args_json["sequence motifs tar"]["counts"]["tfmodisco.seq_contrib.counts.meme."+encid] = {"file.paths": []}
	
	for formats in ["CWM", "CWM-PFM", "hCWM", "hCWM-PFM", "PFM"]:
		meme_file=odir+name+"/merge_folds_new_may_05_24/counts/"+formats
		if os.path.isfile(meme_file):
			args_json["sequence motifs tar"]["counts"]["tfmodisco.seq_contrib.counts.meme."+encid]["file.paths"].append((meme_file,"tfmodisco.seq_contrib.counts."+formats+".meme."+encid))
		else:
			print(meme_file)
			continue
	args_json["sequence motifs tar"]["profile"]["tfmodisco.seq_contrib.profile.meme."+encid] = {"file.paths": []}

	for formats in ["CWM", "CWM-PFM", "hCWM", "hCWM-PFM", "PFM"]:
		meme_file=odir+name+"/merge_folds_new_may_05_24/profile/"+formats
		if os.path.isfile(meme_file):
			args_json["sequence motifs tar"]["profile"]["tfmodisco.seq_contrib.profile.meme."+encid]["file.paths"].append((meme_file,"tfmodisco.seq_contrib.profile."+formats+".meme."+encid))
		else:
			print(meme_file)
			success=False
			break
	if not success:
		continue
			
	if not os.path.isfile(ooutdir+encode_id[name]+".json"):
		f = open(ooutdir+encode_id[name]+".json", "w")
		json.dump(args_json, f, indent=4)
		f.close()


