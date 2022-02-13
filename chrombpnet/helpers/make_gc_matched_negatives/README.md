Scripts to get  non-peak regions gc-matched with a given foreground (or peak) set.

### Usage

```bash 
bash  run.sh [foreground_bed] [exclude_bed] [inputlen] [output_dir] [genome] [genomewide_gc] [fold] [chrom_sizes]
```

The above script runs two python scripts `get_gc_content.py` and `get_gc_matched_negatives.py` and a `bedtools` operation. Briefly `get_gc_content.py` finds the gc content distribution of the given foreground regions. This script filters the foreground regions if `inputlen` region cannot be formed (if the percentage of peaks filtered is high consider reducing the `inputlen` as your genome might be small). The bedtools operation keeps only those genomewide bins that do not fall in exclude bed as candidate negatives for gc-matching. `get_gc_matched_negatives.py` filters the candidate negatives list to as many negative regions as foreround regions such that they are gc-mactched with foreground.

NOTE: Here I assume the script `run.sh` is executed from the main chrombpnet_paper repo - if you want to run these scripts while in this directory modify the paths in `run.sh`

### Input format:

- foreground_bed: Must be a narrowPeak files, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all files. For all these regions we will find the gc-content fraction and will sample negatives that match this distribution.
- exclude_bed: Must be in bed format, with values minimally for chr, start and end. These are the regions that we do not want in our negatives e.g. blacklist regions, peak regions for ATAC/DNASE bias model training etc.
- inputlen: The input length to consider for foreground gc-fraction calculation.
- output_dir: Directory to store the output files. Make sure the directory already exists. The code generates three output files described below in output format section.
- genome: Reference geneome fasta.
- genomewide_gc: Must be TSV file containing chr, start, end and gc_content. These are genome wide regions binned into `inputlen` regions and with a given stride. This pickle file can be generated using the scripts at `get_genomewide_gc_buckets/run.sh`. If you are using hg38 human reference genome and are (1) using a  `inputlen` of 2114 you can use the pre-generated pickle file saved on the cluster `/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_1057.bed` OR (2) if you are using a `inputlen` of 1000 you can use the pre-generated pickle file saved on the cluster `/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_1057.bed`
- fold: A json file containing a dictionary with test,valid and train keys and values with corresponding chromosomes.
- chrom_sizes: TSV file with chromosome name in first column and size in the second column.

### Output format:

The script outs the following three tab seperated files - 

- foreground.gc.bed: This file has 4 columns, chr, start, end and gc fraction of the foreground bed regions centered at the summit. The start to end span a region of `inputlen`.
- candidate.negatives.bed: This file has only candidate negatives (i.e genome wide negatvies after filtering for the exlude set). It has 4 coloumns chr, start, end and gc fraction. The start to end span a region of `inputlen`. 
- negatives.bed: The final negatives with 4 columns with chr, start, end  and gc fractions of the negatives sampled. The start to end span a region of `inputlen`.
