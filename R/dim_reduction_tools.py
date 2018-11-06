import pandas as pd
import numpy as np
import phenograph
#import uncurl 
import dca
import scanpy.api as sc
from MulticoreTSNE import MulticoreTSNE as TSNE
from umap import UMAP
from phate import PHATE
from sklearn.decomposition import TruncatedSVD

# note that we need to explicitly type numeric values coming from R at the moment
# for whatever reason, a lot of stuff that should be an integer value in Python
# is translated into a float

# on further investigation, using **kwargs may get around this typing problem?

def umap_for_r(r_data_frame, 
                n_neighbors=30, 
                min_dist=0.1, 
                n_components=3, 
                metric='euclidean', 
                init = 'spectral', 
                alpha = 1.0, 
                spread = 1.0, 
                bandwidth = 1.0, 
                random_state = None, 
                angular_rp_forest = False, 
                verbose = True):
    umap_embed = UMAP(n_neighbors = int(n_neighbors), min_dist = float(min_dist), n_components = int(n_components), 
        metric = str(metric), init = str(init), alpha = float(alpha), spread = float(spread), bandwidth = float(bandwidth), 
        random_state = random_state, angular_rp_forest = bool(angular_rp_forest))
    umap_cell_embeddings = umap_embed.fit_transform(r_data_frame)
    return umap_cell_embeddings

def mtsne_for_r(r_data_frame, 
                n_components=3, 
                perplexity=30, 
                verbose=1, 
                n_iter=2000, 
                n_jobs=4, 
                angle = 0.5):
    tsne = TSNE(n_components = int(n_components), perplexity = int(perplexity), verbose = int(verbose), 
                n_iter = int(n_iter), n_jobs = int(n_jobs), angle = float(angle))
    multitsne = tsne.fit_transform(r_data_frame)
    return multitsne

def phate_for_r(r_data_frame, 
                n_components = 3, 
                potential_method='sqrt', 
                mds = 'metric', 
                njobs = 4, 
                verbose = True, 
                t=90):
    phate_operator = PHATE(n_components = n_components, potential_method= potential_method, mds = mds, njobs = njobs, verbose = verbose, t=t)
    phate_embeddings = phate_operator.fit_transform(r_data_frame)
    return phate_embeddings
    
def phenograph_for_r(r_data_frame, **kwargs):
    communities, graph, Q = phenograph.cluster(r_data_frame)
    return communities

def fitsne_for_r(r_data_frame,
                no_dims = 3, 
                perplexity = 30.0, 
                theta = 0.5, 
                rand_seed = -1, 
                max_iter = 1000, 
                stop_lying_iter = 200, 
                fft_not_bh=False, 
                ann_not_vptree=False,
                early_exag_coeff=12.0,
                no_momentum_during_exag=False, 
                start_late_exag_iter=-1,
                late_exag_coeff=-1, 
                n_trees=50, 
                search_k=-1, 
                nterms=3, 
                intervals_per_integer=1, 
                min_num_intervals=50, 
                nthreads=0 ):
    r_data_frame = r_data_frame.copy(order="C")
    fitsne_embed = fitsne.FItSNE(r_data_frame, no_dims = int(no_dims), perplexity = float(perplexity), theta = float(theta), 
        rand_seed = int(rand_seed), max_iter = int(max_iter), stop_lying_iter = int(stop_lying_iter), fft_not_bh = bool(fft_not_bh), 
        ann_not_vptree = bool(ann_not_vptree), early_exag_coeff = float(early_exag_coeff), 
        no_momentum_during_exag = bool(no_momentum_during_exag), start_late_exag_iter = int(start_late_exag_iter),
        late_exag_coeff = int(late_exag_coeff), n_trees = int(n_trees), search_k = int(search_k), nterms = int(nterms), 
        intervals_per_integer = int(intervals_per_integer), min_num_intervals = int(min_num_intervals), nthreads = int(nthreads))
    return fitsne_embed

#def uncurl_for_r(r_data_frame, dim_reduct, **kwargs):
def uncurl_for_r(r_data_frame, num_clusters, **kwargs):
    genes = uncurl.max_variance_genes(r_data_frame, nbins = 5, frac = 0.2)
    data_subset = r_data_frame[genes,:]
    M, W, ll = uncurl.run_state_estimation(data_subset, num_clusters, dist='Poiss', disp=False, reps=6, threads = 8)
    tsvd = TruncatedSVD(50)
    #dr = dim_reduct(n_components = 3, **kwargs)
    mw = M.dot(W)
    mw_log = np.log1p(mw)
    mw_tsvd = tsvd.fit_transform(mw_log.T)
    #mw_dr = dr.fit_transform(mw_tsvd)
    return_tuple = (mw_log.T, mw_tsvd)
    return return_tuple
    
def deep_autoencode(seurat_raw_data, cell_names, gene_names):
    adata = sc.AnnData(seurat_raw_data.transpose(), obs=cell_names, var=gene_names)
    adata.obs_names = r.cell_names
    adata.var_names = r.gene_names
    dca_adata = dca.api.dca(adata, threads = 8)
    conv = adata.X.copy(order='C')
    return_vals = (conv, adata.obs_names.values, adata.var_names.values)
    return(return_vals)
