## Generate chrombpnet predictions on given regions and output to a bigwig file. 

The script `chrombpnet_predict_to_bigwig` can be invoked to generate a bigwig track from chrombpnet predictions on a set of input regions. 
Please provide absolute paths to this script. 

```
usage: chrombpnet_predict_to_bigwig [-h] -bm BIAS_MODEL -cm CHROMBPNET_MODEL
                                    -cmb CHROMBPNET_MODEL_NB -r REGIONS -g
                                    GENOME -c CHROM_SIZES -o OUT_PREFIX
                                    [-b BATCH_SIZE] [-t TQDM]
                                    [-d DEBUG_CHR [DEBUG_CHR ...]]


  -bm BIAS_MODEL, --bias-model BIAS_MODEL
                        Path to bias model h5
  -cm CHROMBPNET_MODEL, --chrombpnet-model CHROMBPNET_MODEL
                        Path to chrombpnet model h5
  -cmb CHROMBPNET_MODEL_NB, --chrombpnet-model-nb CHROMBPNET_MODEL_NB
                        Path to chrombpnet model h5
  -r REGIONS, --regions REGIONS
                        10 column BED file of length = N which matches
                        f['projected_shap']['seq'].shape[0]. The ith region in
                        the BED file corresponds to ith entry in importance
                        matrix. If start=2nd col, summit=10th col, then the
                        input regions are assumed to be for
                        [start+summit-(inputlen/2):start+summit+(inputlen/2)].
                        Should not be piped since it is read twice!
  -g GENOME, --genome GENOME
                        Genome fasta
  -c CHROM_SIZES, --chrom-sizes CHROM_SIZES
                        Chromosome sizes 2 column tab-separated file
  -o OUT_PREFIX, --out-prefix OUT_PREFIX
                        Output bigwig file
  -b BATCH_SIZE, --batch-size BATCH_SIZE
  -t TQDM, --tqdm TQDM  Use tqdm. If yes then you need to have it installed.
  -d DEBUG_CHR [DEBUG_CHR ...], --debug-chr DEBUG_CHR [DEBUG_CHR ...]
                        Run for specific chromosomes only (e.g. chr1 chr2) for
                        debugging
```

