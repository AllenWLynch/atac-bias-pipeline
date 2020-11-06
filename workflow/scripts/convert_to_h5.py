

from anndata import AnnData
import pandas as pd
import argparse
from scipy import sparse
import numpy as np

def get_peak_id(data, keys = ['chr','start','end']):
    return data[keys[0]].astype(str) + '_' + data[keys[1]].astype(str).str.strip() + '_' + data[keys[2]].astype(str).str.strip()

def read_bias_file(filepath):
    colnames = 'chr,start,end,barcode,dup_group,bias1,bias2,gc_content,fragment_len,log_duprate,peak_chr,peak_start,peak_end,corrected_count'.split(',')
    fragments = pd.read_csv(filepath, header = None, sep = '\t')
    fragments.columns = colnames
    fragments['peak_id'] = get_peak_id(fragments, keys = ['peak_chr','peak_start','peak_end'])
    return fragments

def read_stats(filepath, columns):
    stats = pd.read_csv(filepath, sep = '\t', header = None)
    stats.columns = columns
    #stats['samples'] = stats.samples.apply(lambda x : list(map(float, x.strip('[]').split(', ')))).apply(np.array)
    return stats

def read_barcode_stats(filepath):
    return read_stats(filepath, ['barcode','mean_log_duprate','fragment_count','samples'])

def read_peak_stats(filepath):
    peak_stats = read_stats(filepath, ['peak_id','median_log_duprate','fragment_count','samples'])
    peak_stats['rank'] = peak_stats.median_log_duprate.rank(method = 'first',na_option = 'keep')
    return peak_stats

def read_sparse_countmatrix(barcode_file, peaks_file, count_matrix_file, min_peak_proportion = 0):

    barcode_stats = read_barcode_stats(barcode_file)
    barcode_stats.index = barcode_stats.index.astype(str)

    peak_stats = read_peak_stats(peaks_file)
    peak_stats.index = peak_stats.index.astype(str)
    
    counts = sparse.load_npz(count_matrix_file)

    data = AnnData(X = counts, obs = barcode_stats, var = peak_stats)

    data.var['accessible_in_cells'] =  np.array((data.X > 0).sum(axis=0)).reshape(-1)

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type = str, required=True)
    parser.add_argument('--var', type = str, required=True)
    parser.add_argument('--X', type = str, required=True)
    parser.add_argument('-o','--output', type = str, required=True)

    args = parser.parse_args()

    matrix = read_sparse_countmatrix(args.obs, args.var, args.X)

    matrix.write(args.output)
