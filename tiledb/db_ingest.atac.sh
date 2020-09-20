#!/bin/bash
db_ingest --tiledb_metadata tier1.encode.atac.tasks.tsv \
	  --array_name db/atac \
	  --overwrite \
	  --chrom_sizes hg38.chrom.sizes \
	  --attribute_config encode_pipeline \
	  --coord_tile_size 10000 \
	  --task_tile_size 1 \
	  --write_chunk 30000000 \
	  --threads 20 \
	  --max_queue_size 50 \
	  --max_mem_g 200

