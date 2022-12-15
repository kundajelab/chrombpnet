import pyBigWig
import numpy as np

def read_chrom_sizes(fname):
    with open(fname) as f:
        gs = [x.strip().split('\t') for x in f]
    gs = [(x[0], int(x[1])) for x in gs if len(x)==2]
 
    return gs

def get_regions(regions_file, seqlen):
    # regions file is assumed to be centered at summit (2nd + 10th column)
    # it is adjusted to be of length seqlen centered at summit

    assert(seqlen%2==0)

    with open(regions_file) as r:
        regions = [x.strip().split('\t') for x in r]

    regions = [[x[0], int(x[1])+int(x[9])-seqlen//2, int(x[1])+int(x[9])+seqlen//2, int(x[1])+int(x[9])] for x in regions]

    return regions

def write_bigwig(data, regions, gs, bw_out, debug_chr=None, use_tqdm=False, outstats_file=None):
    # regions may overlap but as we go in sorted order, at a given position,
    # we will pick the value from the interval whose summit is closest to 
    # current position
    
    chr_to_idx = {}
    for i,x in enumerate(gs):
        chr_to_idx[x[0]] = i

    bw = pyBigWig.open(bw_out, 'w')
    bw.addHeader(gs)
    
    # regions may not be sorted, so get their sorted order
    order_of_regs = sorted(range(len(regions)), key=lambda x:(chr_to_idx[regions[x][0]], regions[x][1]))

    all_entries = []
    cur_chr = ""
    cur_end = 0

    iterator = range(len(order_of_regs))
    if use_tqdm:
        from tqdm import tqdm
        iterator = tqdm(iterator)

    for itr in iterator:
        # subset to chromosome (debugging)
        if debug_chr and regions[i][0]!=debug_chr:
            continue

        i = order_of_regs[itr]
        i_chr, i_start, i_end, i_mid = regions[i]
    
        if i_chr != cur_chr: 
            cur_chr = i_chr
            cur_end = 0
    
        # bring current end to at least start of current region
        if cur_end < i_start:
            cur_end = i_start
    
        assert(regions[i][2]>=cur_end)
    
        # figure out where to stop for this region, get next region
        # which may partially overlap with this one
        next_end = i_end
    
        if itr+1 != len(order_of_regs):
            n = order_of_regs[itr+1]
            next_chr, next_start, _, next_mid = regions[n]
       
            if next_chr == i_chr and next_start < i_end:
                # if next region overlaps with this, end between their midpoints
                next_end = (i_mid+next_mid)//2
       
        vals = data[i][cur_end - i_start:next_end - i_start]

        bw.addEntries([i_chr]*(next_end-cur_end), 
                       list(range(cur_end,next_end)), 
                       ends = list(range(cur_end+1, next_end+1)), 
                       values=[float(x) for x in vals])
    
        all_entries.append(vals)
        
        cur_end = next_end

    bw.close()

    all_entries = np.hstack(all_entries)
    if outstats_file != None:
        with open(outstats_file, 'w') as f:
            f.write("Min\t{:.6f}\n".format(np.min(all_entries)))
            f.write(".1%\t{:.6f}\n".format(np.quantile(all_entries, 0.001)))
            f.write("1%\t{:.6f}\n".format(np.quantile(all_entries, 0.01)))
            f.write("50%\t{:.6f}\n".format(np.quantile(all_entries, 0.5)))
            f.write("99%\t{:.6f}\n".format(np.quantile(all_entries, 0.99)))
            f.write("99.9%\t{:.6f}\n".format(np.quantile(all_entries, 0.999)))
            f.write("99.95%\t{:.6f}\n".format(np.quantile(all_entries, 0.9995)))
            f.write("99.99%\t{:.6f}\n".format(np.quantile(all_entries, 0.9999)))
            f.write("Max\t{:.6f}\n".format(np.max(all_entries)))
