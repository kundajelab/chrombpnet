import json

cell_line="H1ESC"

job_meta = open("hub_template.json").read()

job_new_meta = job_meta.replace("$experiment", cell_line)


#f = open("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/ATAC/"+cell_line+"/hub.json","w")
f = open("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/"+cell_line+"/hub.json","w")

f.write(job_new_meta)
f.close()
