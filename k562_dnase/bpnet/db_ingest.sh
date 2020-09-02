#!/bin/bash
db_ingest --tiledb_metadata K562.tasks.nodup.tsv \
	  --tiledb_group k562_nodup_dnase \
	  --overwrite \
	  --chrom_sizes hg38.chrom.sizes \
	  --attribute_config encode_pipeline \
	  --tile_size 10000 \
	  --write_chunk 10000000 \
	  --chrom_threads 20
