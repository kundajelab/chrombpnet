#!/bin/bash


db_ingest --tiledb_metadata microglia.tsv \
	  --array_name db \
	  --overwrite \
	  --chrom_sizes hg38.chrom.sizes \
	  --attribute_config_file microglia.attribs.txt \
	  --coord_tile_size 10000 \
	  --task_tile_size 1 \
	  --write_chunk 30000000 \
	  --threads 40 \
	  --max_queue_size 50 \
	  --max_mem_g 200

