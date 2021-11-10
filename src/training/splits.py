
chroms = ['chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8', 'chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
splits=dict()

splits[0]={'test':['chr1'],
                'valid':['chr10','chr8'],
                'train':['chr2','chr3','chr4','chr5','chr6','chr7','chr9','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']}

splits[1]={'test':['chr19','chr2'],
                'valid':['chr1'],
                'train':['chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr20','chr21','chr22','chrX','chrY']}

splits[2]={'test':['chr3','chr20'],
                'valid':['chr19','chr2'],
                'train':['chr1','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr21','chr22','chrX','chrY']}

splits[3]={'test':['chr13','chr6','chr22'],
                'valid':['chr3','chr20'],
                'train':['chr1','chr2','chr4','chr5','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr16','chr17','chr18','chr19','chr21','chrX','chrY']}

splits[4]={'test':['chr5','chr16','chrY'],
                'valid':['chr13','chr6','chr22'],
                'train':['chr1','chr2','chr3','chr4','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr17','chr18','chr19','chr20','chr21','chrX']}

splits[5]={'test':['chr4','chr15','chr21'],
                'valid':['chr5','chr16','chrY'],
                'train':['chr1','chr2','chr3','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr17','chr18','chr19','chr20','chr22','chrX']}
    
def get_bed_regions_for_fold_split(bed_regions,fold,split):
    chroms_to_keep=splits[int(fold)][split]
    bed_regions_to_keep=bed_regions[bed_regions["chr"].isin(chroms_to_keep)]
    print("got split:"+str(split)+" for fold:"+str(fold) +" for bed regions:"+str(bed_regions_to_keep.shape))
    return bed_regions_to_keep, chroms_to_keep

    
    
