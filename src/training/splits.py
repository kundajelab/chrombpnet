#generates hg19_splits for cross-validation
#hg19 (i.e. /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes). Any chromosome from this chrom.sizes file that is not in the test/validation split is assumed to be in the training split (only considering chroms 1 - 22, X, Y
hg19_chroms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
hg19_splits=dict()

hg19_splits[0]={'test':['chr1'],
                'valid':['chr10','chr8'],
                'train':['chr2','chr3','chr4','chr5','chr6','chr7','chr9','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']}

hg19_splits[1]={'test':['chr19','chr2'],
                'valid':['chr1'],
                'train':['chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr20','chr21','chr22','chrX','chrY']}

hg19_splits[2]={'test':['chr3','chr20'],
                'valid':['chr19','chr2'],
                'train':['chr1','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr21','chr22','chrX','chrY']}

hg19_splits[3]={'test':['chr13','chr6','chr22'],
                'valid':['chr3','chr20'],
                'train':['chr1','chr2','chr4','chr5','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr16','chr17','chr18','chr19','chr21','chrX','chrY']}

hg19_splits[4]={'test':['chr5','chr16','chrY'],
                'valid':['chr13','chr6','chr22'],
                'train':['chr1','chr2','chr3','chr4','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr17','chr18','chr19','chr20','chr21','chrX']}

hg19_splits[5]={'test':['chr4','chr15','chr21'],
                'valid':['chr5','chr16','chrY'],
                'train':['chr1','chr2','chr3','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr17','chr18','chr19','chr20','chr22','chrX']}

hg19_splits[6]={'test':['chr7','chr18','chr14'],
                'valid':['chr4','chr15','chr21'],
                'train':['chr1','chr2','chr3','chr5','chr6','chr8','chr9','chr10','chr11','chr12','chr13','chr16','chr17','chr19','chr20','chr22','chrX','chrY']}

hg19_splits[7]={'test':['chr11','chr17','chrX'],
                'valid':['chr7','chr18','chr14'],
                'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr8','chr9','chr10','chr12','chr13','chr15','chr16','chr19','chr20','chr21','chr22','chrY']}

hg19_splits[8]={'test':['chr12','chr9'],
                'valid':['chr11','chr17','chrX'],
                'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr10','chr13','chr14','chr15','chr16','chr18','chr19','chr20','chr21','chr22','chrY']}

hg19_splits[9]={'test':['chr10','chr8'],
                'valid':['chr12','chr9'],
                'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr11','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']}


#Note: the splits for hg19 and hg38 are the same, as are the chromosomes used for training models. 
splits=dict()
splits['hg19']=hg19_splits 
splits['hg38']=hg19_splits
chroms=dict()
chroms['hg19']=hg19_chroms
chroms['hg38']=hg19_chroms


mm10_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
mm10_splits = dict()

mm10_splits[0] = {'test':['chr1'],
                  'valid':['chr10','chr8'],
                  'train':['chr2','chr3','chr4','chr5','chr6','chr7','chr9','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']}

mm10_splits[1] = {'test':['chr19','chr2'],
                  'valid':['chr1'],
                  'train':['chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chrX','chrY']}

mm10_splits[2] = {'test':['chr3'],
                  'valid':['chr19','chr2'],
                  'train':['chr1','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chrX','chrY']}

mm10_splits[3] = {'test':['chr13','chr6'],
                  'valid':['chr3'],
                  'train':['chr1','chr2','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr17','chr18','chr19','chrX','chrY']}

mm10_splits[4] = {'test':['chr5','chr16','chrY'],
                  'valid':['chr13','chr6'],
                  'train':['chr1','chr2','chr3','chr4','chr7','chr8','chr9','chr10','chr11','chr12','chr14','chr15','chr17','chr18','chr19','chrX']}

mm10_splits[5] = {'test':['chr4','chr15'],
                  'valid':['chr5','chr16','chrY'],
                  'train':['chr1','chr2','chr3','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr17','chr18','chr19','chrX']}

mm10_splits[6] = {'test':['chr7','chr18','chr14'],
                  'valid':['chr4','chr15'],
                  'train':['chr1','chr2','chr3','chr5','chr6','chr8','chr9','chr10','chr11','chr12','chr13','chr16','chr17','chr19','chrX','chrY']}

mm10_splits[7] = {'test':['chr11','chr17','chrX'],
                  'valid':['chr7','chr18','chr14'],
                  'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr8','chr9','chr10','chr12','chr13','chr15','chr16','chr19','chrY']}

mm10_splits[8] = {'test':['chr12','chr9'],
                  'valid':['chr11','chr17','chrX'],
                  'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr10','chr13','chr14','chr15','chr16','chr18','chr19','chrY']}

mm10_splits[9] = {'test':['chr10','chr8'],
                  'valid':['chr12','chr9'],
                  'train':['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr11','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']}

#mm9 & mm10 splits are the same
splits['mm10']=mm10_splits
splits['mm9']=mm10_splits
chroms['mm10']=mm10_splits
chroms['mm9']=mm10_splits


def get_chroms(args,split):
    if split=="train":
        if args.train_chroms is not None:
            return args.train_chroms
        else:
            assert args.genome is not None
            assert args.fold is not None
            return splits[args.genome][args.fold]['train']
    elif split=='valid':
        if args.validation_chroms is not None:
            return args.validation_chroms
        else:
            assert args.genome is not None
            assert args.fold is not None
            return splits[args.genome][args.fold]['valid']
    elif split=="test":
        if args.predict_chroms is not None:
            return args.predict_chroms
        else:
            assert args.genome is not None
            assert args.fold is not None
            return splits[args.genome][args.fold]['test']
    else:
        raise Exception("invalid split specified, must be one of train,valid,test; you gave:"+str(split))
    

def get_bed_regions_for_fold_split(bed_regions,genome,fold,split):
    chroms_to_keep=splits[genome][int(fold)][split]
    bed_regions_to_keep=bed_regions[bed_regions[0].isin(chroms_to_keep)]
    print("got split:"+str(split)+" for fold:"+str(fold) +" for bed regions:"+str(bed_regions_to_keep.shape))
    return bed_regions_to_keep

    
    
