plink --bfile /oak/stanford/groups/akundaje/refs/1000genomes/hg38_plink_bed/euro.1000G.phase3.hg38.merged \
      --r2 \
      --ld-window-r2 0.8 \
      --allow-extra-chr \
      --ld-snp-list chr1.histoneqtl.rsid.txt
