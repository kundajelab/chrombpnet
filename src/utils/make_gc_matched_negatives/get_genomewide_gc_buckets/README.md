Scripts to find the gc content of genome wide bins.

### Usage

```bash 
python get_genomewide_gc_bins.py -g [reference_fasta]  -c [chromosome_sizes_for_reference] -o [output_path_prefix]  -f [bin_flank_size] --s [stride]
```

The above scripts bins the genome into 2*`flank_size` length regions. The bin intervals are spaced at `stride` length.

### Input format:

- reference_fasta: Reference geneome fasta.
- chromosome_sizes_for_reference: Tab seperated file with two columns. First column is the chromosome and second column is the chromsome length.
- output_path_prefix: Path and file name to store the output. See output fotmat section below.
- bin_flank_size: An integer value for the flank size surrounding the center of a bin for gc-fraction calculation.
- stride: An integer stride value to space the bins.

### Output format:

The script outputs a tab-seperated file containing chr, start, end and gc_content. These are genome wide regions binned into 2*`flank_size` regions and with a given `stride`. 
For example look at the following pre-generated file at `/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_flank_size_1057.bed` generated with `flank_size` 1057 and `stride` 50 using script `/oak/stanford/groups/akundaje/anusri/refs/get_genomewide_gc_bins.sh`.
