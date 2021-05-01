#GM12878
zcat /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz | bedtools sort -i - | bedtools merge -i - > GM12878.overlap.ATAC.DNASE.merged.bed

#H1HESC
zcat /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/58fb3f13-be45-45de-8a39-d0bfbeaf86c5/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz | bedtools sort -i - | bedtools  merge -i - > H1HESC.overlap.ATAC.DNASE.merged.bed 

#K562
bedtools sort -i /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz  | bedtools merge -i - > K562.overlap.DNASE.bed

#HEPG2
zcat /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/ce805260-55f8-43c8-b2a1-a232b4a0e369/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz | bedtools sort -i - | bedtools merge -i - > HEPG2.overlap.ATAC.DNASE.merged.bed 

#IMR90
zcat /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/d754a34e-bc9f-4270-8020-bc37e8d195ba/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/277549db-c2d8-49d3-ace0-81ad5d4088fb/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz | bedtools sort -i - | bedtools merge -i - > IMR90.overlap.ATAC.DNASE.merged.bed
