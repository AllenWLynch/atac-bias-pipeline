
import scanpy as sc
from anndata import AnnData
import pandas as pd
import numpy as np
from scipy import sparse
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize
import seaborn as sns
import subprocess
from functools import partial
from io import StringIO

#__DATA READING__
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
    peak_stats = read_peak_stats(peaks_file)
    
    counts = sparse.load_npz(count_matrix_file)

    data = AnnData(X = counts, obs = barcode_stats, var = peak_stats)

    data.var['accessible_in_cells'] =  np.array((data.X > 0).sum(axis=0)).reshape(-1)

    return data

#__PLOT GENERATION__
def benchmark_fragment_model(fragments, dupcounts, peak_stats):

    fragments = fragments.merge(peak_stats[['peak_id','accessible_in_cells']], on = 'peak_id', how = 'inner')\
        .merge(dupcounts, left_on = 'dup_group', right_on = 'group')
    fragments['true_log_duprate'] = np.log(fragments.duplicate_counts / fragments.accessible_in_cells)

    return fragments[['log_duprate', 'true_log_duprate']]

def sample_ranks(andata_chunk, num_samples = 10000):

    cluster_peak_counts = np.array(andata_chunk.X.sum(axis = 0)).reshape(-1)
    cluster_peak_probs = cluster_peak_counts / cluster_peak_counts.sum()

    peak_samples = np.random.choice(len(cluster_peak_probs), p = cluster_peak_probs, 
        size = (min(len(cluster_peak_probs), num_samples),))
    
    sampled_ranks = andata_chunk.var.iloc[peak_samples]['rank'].values.astype(int)

    return sampled_ranks

def get_bias_peak_enrichment(andata, cluster_key, num_samples = 10000):

    ranks_dict = {}
    for cluster_id in np.unique(andata.obs[cluster_key]):
        ranks_dict[cluster_id] = sample_ranks(
            andata[andata.obs[cluster_key] == cluster_id, ~np.isnan(andata.var['rank'])], 
            num_samples=num_samples
        )

    ranks_df = pd.DataFrame(ranks_dict)
    ranks_df = pd.melt(ranks_df, var_name = 'cluster_id', value_name = 'rank')
    return ranks_df

def get_fragment_distribution_by_peak_and_cluster(andata, peaks_str, cluster_key):

    selected_fragments = read_bias_file(StringIO(peaks_str))

    selected_fragments = selected_fragments.merge(andata.obs[['barcode',cluster_key]], on = 'barcode', how = 'right')\
        .merge(andata.var[['peak_id','rank']], on = 'peak_id', how = 'right')

    stratified_subset = selected_fragments.groupby(['peak_id',cluster_key], group_keys=False)\
        .apply(lambda x : x.sample(min(len(x), 150)))

    return stratified_subset

#__SCANPY ADDONS___
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

def _binary_search(func, acceptible_range, low, high, iter_num = 0, max_iter = 10):

    midpoint = low + (high - low)/2

    output = func(midpoint)
    
    if (output >= acceptible_range[0] and output < acceptible_range[1]) or iter_num >= max_iter:
        return midpoint

    if output < acceptible_range[0]:
        return _binary_search(func, acceptible_range, midpoint, high, iter_num + 1, max_iter=max_iter)
    else:
        return _binary_search(func, acceptible_range, low, midpoint, iter_num+1, max_iter=max_iter)        
        

def get_num_clusters(andata, resolution):
    sc.tl.leiden(andata, resolution = 10**resolution)
    num_clusters = andata.obs.leiden.nunique()
    print('Resolution: {}, # Clusters: {}'.format(str(10**resolution), str(num_clusters)))
    return num_clusters
    
def process_counts(andata, resolution = 0.5, tune = False, num_clusters = (7,10)):
    
    lsi_model = TruncatedSVD(n_components=50)
    X_LSI = lsi_model.fit_transform(
        TfidfTransformer().fit_transform(andata.X)
    )

    andata.obsm['X_lsi'] = X_LSI
    
    sc.pp.neighbors(andata, use_rep = 'X_lsi')

    if tune:
        cluster_func = partial(get_num_clusters, andata)
        resolution = _binary_search(cluster_func, num_clusters, -2, 0)
    else:
        sc.tl.leiden(andata, resolution = resolution)

    sc.tl.umap(andata)
    
    return lsi_model, resolution

def compute_neighborhood_change(*,control, treatment, cluster_key='leiden', mask_cluster = True):
    
    a1_distances = cosine_similarity(control.obsm['X_lsi'])
    a2_distances = cosine_similarity(treatment.obsm['X_lsi'])
    
    if mask_cluster:
        clusters = control.obs[cluster_key].values
        cluster_square = np.tile(clusters, (len(clusters), 1))

        same_cluster_mask = cluster_square == cluster_square.T
    else:
        same_cluster_mask = np.full(a1_distances.shape, True)
    
    a1_distances = normalize(np.where(same_cluster_mask, 0.0, a1_distances), 'l2')
    a2_distances = normalize(np.where(same_cluster_mask, 0.0, a2_distances), 'l2')
    
    neighborhood_change = 1 - np.sum(np.multiply(a1_distances, a2_distances), axis = 1)
    
    control.obs['neighbor_change'] = neighborhood_change