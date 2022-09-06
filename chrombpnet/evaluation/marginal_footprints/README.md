# Scripts to do marginal footprinting using ChromBPNet

The scripts in this folder provide marginal footprints for a given motif and background regions.  Briefly the algorithm can be desicribed in the following three steps -  (1) Insert a given motif's center in the center of background regions to make synthetic sequences. (2) Find profile probability predictions for the synthetic sequences created for a given motif. (3) Find profile probability predictions for the reverse complement of the given synthetic sequences and then complement the predicionts. (4) Average the prediction in (2) and (3) to get footprinting for a given synthetic sequence. (3) Average the footprints across all the synthetic sequences to get marginal footprint for a given motif.

## Usage

```
chrombpnet_marginal_footprints -g [genome_fasta] -r [regions] -chr [test_chr] -m [model_h5] -bs [batch_size] -o [output_prefix] -pwm_f [motifs_to_pwm_file] --ylim (min-y, max-y)
```

## Example Usage

```
chrombpnet_marginal_footprints -g /path/to/genome_fasta -r /path/to/bed_file -chr chr1 -m /path/to/model.h5 -o /path/to/output_dir/outputprefix -pwm_f /path/to/motif_to_pwm.tsv
```

## Input Format

- genome_fasta: Reference geneome fasta.
- regions: Must be a narrowPeak files, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all files. 
- test_chr: Chromosome name to use. The name mentioned here should be present in column 1 in the `regions` bed file.  We only use regions corresponding to this chromosome.
- model_h5: Model in hdf5 format.
- batch_size: Batch size to use for model predictions.
- output_prefix: Output prefix path and name to use. See Output format section below to understand how this is used.
--motifs_to_pwm_file:  Path to a TSV file containing motifs in first column (e.g. `Tn5`) and motif string (e.g. `GCACAGTACAGAGCTG`) to use for footprinting in second column. A default file is provided in the data folder (https://github.com/kundajelab/chrombpnet/tree/master/data/motif_to_pwm.tsv)
--ylim: optional argument to specify the lower and upper y-axis values for the generated plot. This is a tuple of float/int values (lower y-lim bound for plotting, upper-ylim bound for plotting). Default is to not use this argument and auto-calculate the y-axis range.

Example  https://github.com/kundajelab/chrombpnet/tree/master/data/motif_to_pwm.tsv can be used as a default and contains the following:
```
tn5_1    GCACAGTACAGAGCTG
tn5_2    GTGCACAGTTCTAGAGTGTGCAG
tn5_3    CCTCTACACTGTGCAGAA
tn5_4    GCACAGTTCTAGACTGTGCAG
tn5_5    CTGCACAGTGTAGAGTTGTGC
dnase_1    TTTACAAGTCCA
dnase_2    TGTACTTACGAA
NRF1    GCGCATGCGC
AP1    CGATATGACTCATCCC
CTCF    TTGGCCACTAGGGGGCGCTAT
ETS    CCGAAAGCGGAAGTGAGAC
SP1    AAGGGGGCGGGGCCTAA
RUNX    CCCTAACCACAGCCC
NFKB    GCAAGGGAAATTCCCCAGG
GATA+TAL    GGCTGGGGGGGGCAGATAAGGCC
TAL    GGCTGGG
NFYB    CCAGCCAATCAGAGC
GABPA    GAAACCGGAAGTGGCC
BACH1+MAFK    AACTGCTGAGTCATCCCG
NRF1    CCCCGCGCATGCGCAGTGC
HNF4G    CCGTTGGACTTTGGACCCTG
```


## Output Format

The following two files are created using the `output_prefix` as prefix for the output  Note prefix xan include a directory oath and name.

- `output_prefix`.footprints.h5: A h5 formatted file containing a dictionary with `motifs` as key names and average `model_h5` predictions as values.
- `output_prefix`.`motif_name`.footprints.png: A set of png images each with the naming convention as follows - `output_prefix`.`motif_name`.footprints.png. `motif_name` is a value from the list of `motifs` input. And the image carries the center 200bp marginal footprint for that motif.

