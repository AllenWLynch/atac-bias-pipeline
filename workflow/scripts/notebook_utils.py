
import scanpy as sc
from anndata import AnnData
import pandas as pd
import numpy as np
from scipy import sparse
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
import seaborn as sns

def get_peak_id(data, keys = ['chr','start','end']):
    return data[keys[0]].astype(str) + '_' + data[keys[1]].astype(str).str.strip() + '_' + data[keys[2]].astype(str).str.strip()

def read_bias_file(filepath):
    colnames = 'chr,start,end,barcode,dup_group,bias1,bias2,gc_content,fragment_len,log_duprate,peak_chr,peak_start,peak_end,corrected_count'.split(',')
    fragments = pd.read_csv(filepath, header = None, sep = '\t')
    fragments.columns = colnames
    fragments['peak_id'] = get_peak_id(fragments, keys = ['peak_chr','peak_start','peak_end'])
    return fragments

def benchmark_fragment_model(fragments):

    dup_group_counts = fragments['dup_group'].value_counts()
    dup_group_counts.name = 'duplicate_counts'
    dup_group_counts[-1] = 1

    peak_counts = fragments.groupby('peak_id')['barcode'].nunique()
    peak_counts.name = 'peak_accessible_in_cells'

    fragments = fragments.sample(5000)
    fragments = fragments.merge(peak_counts, on = 'peak_id', how = 'inner')\
        .merge(dup_group_counts.reset_index(), left_on = 'dup_group', right_on = 'index')
    fragments['true_log_duprate'] = np.log(fragments.duplicate_counts / fragments.peak_accessible_in_cells)

    return fragments[['log_duprate', 'true_log_duprate']]

def aggregate_cell_stats(fragments):

    groupby_object = fragments.groupby('barcode')

    barcode_bias = groupby_object.agg({'log_duprate' : np.nanmean, 'bias1' : 'count'})
    barcode_bias.columns = ['mean_log_duprate', 'fragment_count']
    barcode_list = groupby_object['log_duprate'].apply(list)
    barcode_list.name = 'all_log_duprates'

    return barcode_bias.join(barcode_list)

def rank_peaks(fragments):

    peak_ranks = fragments.groupby('peak_id')['log_duprate'].median().sort_values().reset_index()\
        .reset_index().rename(columns = {'index' : 'rank'})

    peak_ranks = peak_ranks.merge(fragments.peak_id.value_counts().reset_index()\
        .rename(columns = {'peak_id' : 'fragment_count','index' : 'peak_id'}), on = 'peak_id')

    return peak_ranks

def get_bias_peak_enrichment(fragments, clusters, peak_ranks, num_samples = 10000):

    clusters.name = 'cluster'

    cluster_ranks = fragments[['barcode','peak_id']].merge(peak_ranks[['peak_id', 'rank']], on = 'peak_id')\
        .merge(clusters.reset_index().rename(columns = {'index' : 'barcode'}), on = 'barcode', how = 'left').groupby('cluster')['rank']\
        .apply(lambda x : np.random.choice(x.values.reshape(-1), size = min(num_samples, len(x)), replace = False)).explode()
    cluster_ranks.columns = ['cluster','rank']

    return cluster_ranks

def get_fragment_distribution_by_peak_and_cluster(fragments, clusters, peak_ranks):

    unbiased_rank, biased_rank = np.quantile(peak_ranks['rank'].values, [0.1,0.9])
    peaks_of_interest = pd.concat([
        peak_ranks[peak_ranks['rank'] >= biased_rank].sample(10),
        peak_ranks[peak_ranks['rank'] <= unbiased_rank].sample(10)
    ])

    clusters.name = 'cluster'

    fragments = peaks_of_interest.merge(fragments, on = 'peak_id', how = 'left')\
        .merge(clusters.reset_index().rename(columns = {'index' : 'barcode'}), on = 'barcode', how = 'right')

    stratified_subset = fragments.groupby(['peak_id','cluster'], group_keys=False)\
        .apply(lambda x : x.sample(min(len(x), 150))).sort_values('rank')

    return stratified_subset

def read_sparse_countmatrix(barcode_file, peaks_file, count_matrix_file, min_peak_proportion = 0):

    with open(barcode_file, 'r') as f:
        barcodes = [b.strip() for b in f]

    with open(peaks_file, 'r') as f:
        peaks = [p.strip().replace('\t','_') for p in f]

    barcodes = pd.DataFrame(index = barcodes)
    peaks = pd.DataFrame(index = peaks)

    counts = sparse.load_npz(count_matrix_file)

    data = AnnData(X = counts, obs = barcodes, var = peaks)
    data.obs['unique_peaks'] = (data.X > 0).sum(axis = 1)
    data = data[data.obs.unique_peaks >= min_peak_proportion * len(peaks)]

    return data

def get_differential_peaks(andata, top_n = 200):

    sc.tl.rank_genes_groups(andata, 'leiden', method='wilcoxon')

    return pd.DataFrame(andata.uns['rank_genes_groups']['names']).head(200)

def standardize(x):
    return np.clip((x - x.mean()) / x.std(), -3, 3)

def plot_umap(andata,*,color_key, key = 'X_umap', quantitative=True, ax = None, center=True, legend = False,
             hue_order = None, hue_norm = None):
    data = andata.obsm['X_umap']
    if type(color_key) == str:
        color_var = andata.obs[color_key].values
    else:
        color_var = color_key
        
    if not quantitative:
        ax = sns.scatterplot(x = data[:,0], y = data[:,1], size = 0.5, 
                       hue = color_var.astype('str'), 
                       hue_order = hue_order, legend = legend, ax = ax)
    else:
        ax = sns.scatterplot(x = data[:,0], y = data[:,1], hue = standardize(color_var) if center else color_var,
                   palette = "vlag", size = 0.5, legend = False, ax = ax, hue_norm = hue_norm)
    ax.set(xticklabels = [], xticks = [], yticks = [])
    sns.despine()
    return ax
    
def process_counts(andata):
    
    lsi_model = TruncatedSVD(n_components=50)
    X_LSI = lsi_model.fit_transform(
        TfidfTransformer().fit_transform(andata.X)
    )

    andata.obsm['X_lsi'] = X_LSI
    
    sc.pp.neighbors(andata, use_rep = 'X_lsi')

    sc.tl.leiden(andata, resolution = 0.5)

    sc.tl.umap(andata)
    
    return lsi_model