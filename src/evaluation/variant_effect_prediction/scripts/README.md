
# Scripts to do variant effect prediction using ChromBPNet

The following scripts in this folder are used to do variant effect prediction using chrombpnet models. We tested these models on bQTLS (Tehranchi et al. 2016) and dsQTLS (Degner et al 2012). We benchmarked the performance of dsQTLs with deltaSVM (Lee et al. 2015). Scripts to reproduce these results will be provided soon.

## Usage

```
python snp_scoring.py -i [snp_data_tsv] -g [genome_fasta] -m [model_hdf5] -o [output_dir] -bs [batch_size]
```

The following assumptions are made with this script - make changes accordingly if the assumptions dont hold.

- The script is designed to work only with SNPs - so make sure reference and alternate allele are just one character. 
- The script returns the following three effect scores - 
    - log counts difference in alternate allele and reference allele predictions (`log_counts_diff`)
    - Sum of absolute difference in log probabilites per base (`log_probs_diff_abs_sum`)
    - Jensenshannon distance between alternate allele profile probability predictions and reference allele profile probability predictions (`probs_jsd_diff`).
- Carefully read the input format of `snp_data_tsv` below.
- The input sequence length (`inputlen`) is inferred from the model. In chrombpnet models the `inputlen` used is an even number to ensure symmtery. So here we insert the allele at `inputlen`//2 locus. Which means that the sequence left of allele is 1bp longer than the sequence right of allelle.
- If the reference/alternate allele are at the edge of chromosome - preventing us from generating an `inputlen` sequence - we will ignore this SNP and print the locus information of the ignored SNP. This might result in the final output snp predicted being smaller than the given input snps.
- The reference genome loader used is `pyfaidx` which is 0-based and the allele position provided in `snp_data_tsv` is expected to be 0-based too. 
- If unsure wheather you are inserting the allele correctly - verify that the sequence generate for a given allele is what is expected to be found by referring to https://www.ncbi.nlm.nih.gov/snp/.

## Example Usage

```
python snp_scoring.py -i /mnt/lab_data2/anusri/variant_effect_prediction_example/subsample_test.csv -g /mnt/data/male.hg19.fa -m /path/to/model.hdf5 -o /path/to/store/output -bs 64
```

## Input Format

- snp_data_tsv: A TSV file with the following 5 columns -  chr, position (0-based) to insert the allele, reference allele, alternate alllele, meta-data. You can leave the meta-data empty too. 
    - Meta-data column can be used to provide the following information such as pvalue significance of the snp, observed effect scores etc. This information will be added as column information to the output `variant_scores.tsv` (see below). For example look at the input tsv file that I use here - `/mnt/lab_data2/nusri/variant_effect_prediction_example/subsample_test.csv` where the meta-data provides information such as significance and observed effect scores as comma-seperated value.
    - Make sure that there are no duplicates in this file - that is the following 4 columns -  chr, position (0-based) to insert the allele, reference allele, alternate alllele, - are not repeated with different meta-data values.
- genome_fasta: Reference geneome fasta. 
- model_hdf5: Model in hdf5 format.
- output_dir: Directory to store the output files. Make sure the directory already exists. The code generates two output files described below in output format section.
- batch_size: Batch size to use for model predictions.

## Output Format

variant_scores.tsv: A TSV file with 8 columns - chr, position (0-based) where the allele was inserted, reference allele, alternate alllele, `log_counts_diff`,  `log_probs_diff_abs_sum`, `probs_jsd_diff` and meta-data. The meta-data if presented in `snp_data_tsv` will be copied to this file, otherwise it will have empty column with no values per row.
predictions_at_snp.pkl: A pickle file containing a dictionary with the following keys - `rsids`, `ref_logcount_preds`, `alt_logcount_preds`, `ref_prob_preds`, `alt_prob_preds`. 