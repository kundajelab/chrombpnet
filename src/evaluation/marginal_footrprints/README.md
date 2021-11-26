# Scripts to do marginal footprinting using ChromBPNet

The scripts in this folder provide marginal footprints for a given motif and background regions.  Briefly the algorithm can be desicribed in the following three steps -  (1) Insert a given motif in the center of background regions to make synthetic sequences. (2) Find profile probability predictions for the synthetic sequences created for a given motif. (3) Find profile probability predictions for the reverse complement of the given synthetic sequences and then complement the predicionts. (4) Average the prediction in (2) and (3) to get footprinting for a given synthetic sequence. (3) Average the footprints across all the synthetic sequences to get marginal footprint for a given motif.

## Usage

```
python marginal_footprinting.py -g [genome_fasta] -r [regions] -chr [test_chr] -m [model_h5] -bs [batch_size] -o [output_prefix] -pwm_f [motifs_to_pwm_file] -mo [motifs]
```

## Example Usage

```
python marginal_footprinting.py -g /path/to/genome_fasta -r /path/to/bed_file -chr chr1 -m /path/to/model.h5 -o /path/to/output_dir/outputprefix -pwm_f motif_to_pwm.tsv -mo CTCF Tn5
```

## Input Format

- genome_fasta: Reference geneome fasta.
- regions: Must be a narrowPeak files, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all files. 
- test_chr: Chromosome name to use. The name mentioned here should be present in column 1 in the `regions` bed file.  We only use regions corresponding to this chromosome.
- model_h5: Model in hdf5 format.
- batch_size: Batch size to use for model predictions.
- output_prefix: Output prefix path and name to use. See Output format section below to understand how this is used.
- motifs_to_pwm_file: Path to a TSV file containing motifs in first column (e.g. `Tn5`) and motif string (e.g. `GCACAGTACAGAGCTG`) to use for footprinting in second column. 
- motifs: The motifs to filter from `motifs_to_pwm_file` for footprinting. We will find footprints only for the motifs mentioned here. Make sure that the motif names mentioned here are present in the column 1 of `motifs_to_pwm_file`. Provide more than one motif - seperate them by space.


## Output Format

The following two files are created using the `output_prefix` as prefix for the output  Note prefix xan include a directory oath and name.

- `output_prefix`.footprints.h5: A h5 formatted file containing a dictionary with `motifs` as key names and average `model_h5` predictions as values.
- `output_prefix`.`motif_name`.footprints.png: A set of png images each with the naming convention as follows - `output_prefix`.`motif_name`.footprints.png. `motif_name` is a value from the list of `motifs` input. And the image carries the center 200bp marginal footprint for that motif.
