## How to run?

```
python hitcalling.py \
        -g [reference genome (fasta file)] \
        -t [tomtom annotation file for motifs (TSV, optional)] \
        -r [genomic regions to run hit caller (bed format - 10 columns)] \
        -cbw [contributions at genomic regions (bigwig)] \
        -bgr  [gc-matched genomic regions to run hit caller for background (bed format - 10 columns)] \
        -bgcbw [contributions at gc-matched genomic regions for background (h5py)] \
        -mo [modisco h5py objects] \
        -o [output dir]  \
        --modisco-motifs-exclude [.txt file for motifs to exclude (for example to exclude `metacluster_0_pattern_0` and  `metacluster_0_pattern_1` for hit call annotations provide "0_0" and "0_1" per line)] \
        --debug
```

## What motifs to exclude using `modisco-motifs-exclude` ?

Exclude any noisy modisco motifs, that have high GC/AT content. 

## Output format

The script creates a bed file `hit_calls.bed` in the output dir. The column names are as follows - 

"chrom", "start", "end", "hit_name", "strand", "pwm_match_score", "mean_contribs", "corr_cwm_with_hit", "p-value", "q-value"
