import pandas as pd
import matplotlib.pyplot as plt

cols=[
        "chrom", "start", "end", "key", "strand", "peak_index",
        "imp_total_signed_score", "imp_total_score", "imp_frac_score",
        "imp_ic_avg_score", "agg_sim", "mod_delta", "mod_precision",
        "mod_percentile", "fann_perclasssum_perc", "fann_perclassavg_perc"
    ]

data = pd.read_csv("results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret/tfhits/tfm_matches.bed", sep="\t", header=None, names=cols)
new_data = data[data["key"]=="0_1"].sort_values(by=['imp_total_score'],ascending=False)
new_data.to_csv("temp_tf.tsv", sep="\t")
print(data[data["key"]=="0_1"]["imp_total_score"].values)
plt.figure()
plt.hist(data[data["key"]=="0_1"]["imp_total_score"].values)
plt.show()
