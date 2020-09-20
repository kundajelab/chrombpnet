#!/bin/bash
db_ingest --tiledb_metadata tier1.encode.chip.tasks.tsv \
	  --array_name db/histone2 \
	  --chrom_sizes hg38.chrom.sizes \
	  --attribute_config_file tier1.encode.chip.attribs.txt \
	  --coord_tile_size 10000 \
	  --task_tile_size 1 \
	  --write_chunk 10000000 \
	  --threads 40 \
	  --max_queue_size 50 \
	  --max_mem_g 200

