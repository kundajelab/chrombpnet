## 1 Generate shifted tagAlign from Granges object

wget http://mitra.stanford.edu/kundaje/projects/chrombpnet/microglia/inputs/Cluster24-fragments.rds

```
    Rscript granges_to_tagAlign.R Cluster24-fragments.rds Cluster24.tagAlign.gz 99
```

## 2. Run the resulting shifted tagAlign file through caper pipeline to generate peak calls

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
