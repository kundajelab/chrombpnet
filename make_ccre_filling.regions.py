import os
import pandas as pd

#experiments=["GM12878_new", "IMR90_new", "HEPG2", "H1ESC_new", "K562"]
experiments=[ "H1ESC_new", "K562"]
encid={"GM12878_new":"ENCSR000EMT",
	"HEPG2": "ENCSR149XIL",
	"IMR90_new": "ENCSR477RTP",
	"K562": "ENCSR000EOT",
	"H1ESC_new": "ENCSR000EMU"}

def get_ccre_regions(experiment, peaks_to_fill, final_regions,  fdr_merged_bed, prep_commands, errors, bed_dir_int, no_ccres_list, log_file):
	url="http://users.wenglab.org/moorej3/Registry-cCREs-WG/V4-Files/"
	ccre_temp_file_name=os.path.join(bed_dir_int,"cre.bed.gz")
	ccre_regions_preformatted=os.path.join(bed_dir_int,"cre.reformat.bed")
	
	celld = pd.read_csv("ccres_cell_type_specfic.csv",sep=",")
	val = celld[celld["Experiment"]==experiment]["cCRE File"].values
	if len(val) == 1:
		command = "wget "+url+val[0]+" "+"-O "+ccre_temp_file_name
		print(command)
		prep_commands.append(command)
		exit_code = os.system(command)
		if exit_code > 0:		
			prep_commands.append(str(exit_code))
			log_lines = "\n".join(prep_commands)
			f = open(log_file, "w")
			f.write(log_lines)
			f.close()
			errors.append([experiment, "error in wget"])
			return False, no_ccres_list, prep_commands, errors

		command = "bedtools intersect -v -wa -a "+ccre_temp_file_name+" -b "+fdr_merged_bed+" | sort | uniq | grep -v \"Low-DNase\" > "+ccre_regions_preformatted
		print(command)
		exit_code = os.system(command)
		prep_commands.append(command)
		if exit_code > 0:		
			prep_commands.append(str(exit_code))
			log_lines = "\n".join(prep_commands)
			f = open(log_file, "w")
			f.write(log_lines)
			f.close()
			errors.append([experiment, "error in bedtools intersect -v ccre"])
			return False, no_ccres_list, prep_commands, errors

		command = "rm "+ccre_temp_file_name
		os.system(command)
		prep_commands.append(command)

		# check if file is non empty
		if os.stat(ccre_regions_preformatted).st_size > 0:
			bed_ccre = pd.read_csv(ccre_regions_preformatted,sep="\t",header=None)
			print(bed_ccre.shape[0])
			if bed_ccre.shape[0] >= 1:
				bed_ccre_new=pd.DataFrame()
				bed_ccre_new[0] = bed_ccre[0]
				bed_ccre_new[1] = bed_ccre[1]
				bed_ccre_new[2] = bed_ccre[2]
				bed_ccre_new[3] = "."
				bed_ccre_new[4] = "."
				bed_ccre_new[5] = "."
				bed_ccre_new[6] = "."
				bed_ccre_new[7] = "."
				bed_ccre_new[8] = "."
				bed_ccre_new[9] = ((bed_ccre[1]+bed_ccre[2])//2)-bed_ccre[1]

#				if os.stat(peaks_to_fill).st_size > 0:
				if peaks_to_fill is not None:

						command = "combine in python "+peaks_to_fill+" "+ccre_regions_preformatted
						prep_commands.append(command)
						old_bed = pd.read_csv(peaks_to_fill,sep="\t",header=None)
						print(old_bed.shape)
						final_bed = pd.concat([bed_ccre_new, old_bed])
				else:
						command = "empty peaks_to_fill use only"+ccre_regions_preformatted
						prep_commands.append(command)
						final_bed = bed_ccre_new
				final_bed.to_csv(final_regions, header=False, index=False, sep="\t", compression='gzip')

				command = "rm "+ccre_regions_preformatted
				os.system(command)
				prep_commands.append(command)

				return True, no_ccres_list, prep_commands, errors
		else:
			#if os.stat(peaks_to_fill).st_size > 0:
			if peaks_to_fill is not None:
				command = "empty "+ccre_regions_preformatted+" use only"+peaks_to_fill
				prep_commands.append(command)
				old_bed = pd.read_csv(peaks_to_fill, sep="\t", header=None)
				final_bed = old_bed
			else:
				command = "both empty "+ccre_regions_preformatted+" and"+peaks_to_fill
				prep_commands.append(command)
				final_bed = pd.DataFrame()
			final_bed.to_csv(final_regions, header=False, index=False, sep="\t", compression='gzip')
			
			command = "rm "+ccre_regions_preformatted
			os.system(command)
			prep_commands.append(command)

			return True, no_ccres_list, prep_commands, errors
	else:
		print("no ccre file")

		no_ccres_list.append(experiment)
		command = "no html for "+ccre_regions_preformatted+" use only"+peaks_to_fill
		prep_commands.append(command)
		#command = "mv "+peaks_to_fill+" "+final_regions
		if os.stat(peaks_to_fill).st_size > 0:
			final_bed = pd.read_csv(peaks_to_fill, sep="\t", header=None)
			
			final_bed.to_csv(final_regions, header=False, index=False, sep="\t", compression='gzip')
		else:
			command = "both empty "+ccre_regions_preformatted+" and"+peaks_to_fill
			prep_commands.append(command)
			final_bed = pd.DataFrame()
			final_bed.to_csv(final_regions, header=False, index=False, sep="\t", compression='gzip')
		return True, no_ccres_list, prep_commands, errors
	return

for cellt in experiments:

	experiment=encid[cellt]
	print(experiment, cellt)

	fdr_merged_bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cellt+"/preds_upload/average_preds/merged.viz.bed.gz"
	peaks_to_fill=None
	final_regions="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cellt+"/ccres/ccres_fill.bed.gz"
	prep_commands=[]
	errors=[]
	bed_dir_int="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cellt+"/ccres/"
	no_ccres_list=[]
	log_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cellt+"/ccres/log.txt"

	get_ccre_regions(experiment, peaks_to_fill, final_regions,  fdr_merged_bed, prep_commands, errors, bed_dir_int, no_ccres_list, log_file)
