## 1 Generate shifted tagAlign from Granges object

wget http://mitra.stanford.edu/kundaje/projects/chrombpnet/microglia/inputs/Cluster24-fragments.rds

```
    Rscript granges_to_tagAlign.R Cluster24-fragments.rds Cluster24.tagAlign.gz 99
```

## 2. Run the resulting shifted tagAlign file through caper pipeline to generate peak calls
Using json with caper: 

```
{
    "atac.pipeline_type": "atac",
    "atac.title": "cluster24_shifted",
    "atac.paired_end": true,
    "atac.genome_tsv": "s3://caper-in/refs/hg38/hg38.tsv",
    "atac.align_cpu": 10,
    "atac.align_mem_factor":0.50,
    "atac.align_disk_factor":20,
    "atac.tas": [
        "s3://caper-in/inbox/Cluster24.tagAlign.sorted.shifted.bed.gz"
    ]
}
```
Outputs are here:
Using version 2.0.3 of pipeline 
```
http://mitra.stanford.edu/kundaje/projects/chrombpnet/microglia/inputs/caper_out/
```

## 3. generate tiledb inputs
```
tiledb/db_ingest_microglia.sh
```


## 4. Train bias model + sequence to counts model. 

```
microglia_atac_bias_filters_128.sh
microglia_atac_bias_filters_500.sh
```

## 5. Run modisco
```
run_modisco_128.sh
run_modisco_500.sh  
run_modisco_bias_128.sh  
run_modisco_bias_500.sh  
trim_modisco.sh  
fetch_tomtom_hits.sh  
aggregate_tomtom_hits.sh  
make_html_modisco_vis.sh  
```
