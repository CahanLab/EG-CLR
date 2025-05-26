import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import matplotlib.pyplot as plt
from muon import atac as ac
from muon import prot as pt
import anndata as ad
import random
random.seed(10)
from scipy import sparse
import pybedtools
import pickle
import anndata as ad
import pySingleCellNet as pySCN # pip install git+https://github.com/pcahan1/PySingleCellNet.
import random
random.seed(10)
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import normalized_mutual_info_score
from sklearn.feature_selection import mutual_info_regression
from scipy import stats

# ---------------------------------------------Utils Functions------------------------------------------------


# limit adata to the shared cells #
def find_common_cell(adata_rna, adata_atac):
    
    # obtain adata cell names
    a = adata_rna.obs_names
    b = adata_atac.obs_names
    
    # Find the intersection of names
    common_cells = b.intersection(a)    
    
    # select cells
    adata_common_RNA = adata_rna[common_cells].copy()
    adata_common_ATAC = adata_atac[common_cells].copy()
    
    print("adata_RNA:",len(adata_common_RNA), ", adata_ATAC:", len(adata_common_ATAC) )
    
    return adata_common_RNA, adata_common_ATAC


# limit to only chr values for atac peaks #
def limit_to_Chr(adata, key = 'Chromosome', char = 'chr'):
    return adata[:,adata.var[key].str.startswith(char,na=False)].copy()


# limit adata to the shared genes 
def find_common_gene(adata_rna, adata_atac, rna_key, atac_key):
       
    # Extract gene names from 'adata'
    rna_gene_names = list(adata_rna.var[rna_key])  
    atac_gene_names = list(adata_atac.var[atac_key])  
        
    # Find the intersection of gene names
    overlap_genes = set(atac_gene_names).intersection(set(rna_gene_names))

    # Subset 'adata' to keep only the overlapping genes
    rna_list = adata_rna.var_names[adata_rna.var[rna_key].isin(overlap_genes)]  
    atac_list = adata_atac.var_names[adata_atac.var[atac_key].isin(overlap_genes)]     
    
    # select genes
    adata_rna_flitered = adata_rna[:,rna_list]
    adata_atac_flitered = adata_atac[:,atac_list]
        
    return adata_rna_flitered, adata_atac_flitered


# limit genes to a file of gene names 
def limit_gene_by_file(file_path, adata, var_key):
    
    # Read gene names from the text file
    with open(file_path, 'r') as file:
        file_gene_names = [line.strip() for line in file]
    
    # Extract gene names from 'adata'
    adata_gene_names = list(adata.var[var_key])  
    
    # Find the intersection of gene names
    overlap_genes = set(file_gene_names).intersection(set(adata_gene_names))

    # Subset 'adata' to keep only the overlapping genes
    adata_TF = adata.var_names[adata.var[var_key].isin(overlap_genes)]  
    adata_filtered = adata[:,adata_TF].copy()
        
    return adata_filtered


# limit genes to a list of gene names 
def limit_gene_by_list(gene_names, adata, var_key):
    
    # Extract gene names from 'adata'
    adata_gene_names = list(adata.var[var_key])  
    
    # Find the intersection of gene names
    overlap_genes = set(gene_names).intersection(set(adata_gene_names))

    # Subset 'adata' to keep only the overlapping genes
    adata_TF = adata.var_names[adata.var[var_key].isin(overlap_genes)]  
    adata_filtered = adata[:,adata_TF]
        
    return adata_filtered


# input adata, save gene and their corresponding values
def df_gene_chromosome(adata):
    
    # select desired to keep col names
    var_df = adata.var
    var_df = var_df[["Chromosome", "gene", "Start", "End", "gene_ids"]]
    
    # drop duplciated genes 
    var_df = var_df.sort_values(by='gene')
    var_df.index = var_df.gene
    var_df = var_df[~var_df.index.duplicated(keep='first')]
    
    return var_df


# write CLR matrixes into pickle file #
def write_matrixes(file_name, matrixes):
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(matrixes, f)


# load CLR matrixes pickle file #
def load_matrxies(file_name):
    with open(file_name, 'rb') as file:
        matrixes = pickle.load(file)
    return matrixes


# ---------------------------------------------Utils Functions------------------------------------------------
