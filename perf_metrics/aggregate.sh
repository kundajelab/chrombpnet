#baseline models 
kerasAC_aggregate_summaries --summary_stats_file baseline.summaries.txt \
			    --out baseline.metrics.txt

#bias corrected models 
kerasAC_aggregate_summaries --summary_stats_file bias.corrected.summaries.txt \
			    --out bias.corrected.metrics.txt 

#pred from bias 
kerasAC_aggregate_summaries --summary_stats_file pred.from.bias.summaries.txt \
			    --out pred.from.bias.metrics.txt 

