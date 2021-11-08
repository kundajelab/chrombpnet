#!/bin/bash

db_ingest --tiledb_metadata $1 \
	  --array_name $2 \
	  --overwrite \
	  --chrom_sizes $3 \
	  --attribute_config_file $4 \
	  --coord_tile_size 10000 \
	  --task_tile_size 1 \
	  --write_chunk 30000000 \
	  --threads 20 \
	  --max_queue_size 50 \
	  --max_mem_g 200

