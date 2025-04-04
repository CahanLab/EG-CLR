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


# write CLR matrixes into pickle file
def write_matrixes(file_name, matrixes):
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(matrixes, f)


# load CLR matrixes pickle file
def load_matrxies(file_name):
    with open(file_name, 'rb') as file:
        matrixes = pickle.load(file)
    return matrixes


# ---------------------------------------------Utils Functions------------------------------------------------





# ---------------------------------------------ATAC-preprocessing Functions------------------------------------------------


# sign chromosome location to gene by taking feature.tsv.gz from filter/raw matrix folder #
def assign_chr(path, adata):
    
    # Load the full features.tsv file including all columns
    features_df = pd.read_csv(path, header=None, sep='\t',compression="gzip")
    features_df.columns = ['gene_ids', 'gene_symbols', 'feature_type', 'Chromosome', 'Start', 'End']
    features_df.set_index("gene_ids", inplace=True)
    
    chromosome_mapping = features_df['Chromosome']
    start_mapping = features_df['Start']
    end_mapping = features_df['End']
    

    # Map the Chromosome data to adata.var using the gene_ids
    adata.var['Chromosome'] = adata.var['gene_ids'].map(chromosome_mapping)
    adata.var['Start'] = adata.var['gene_ids'].map(start_mapping)
    adata.var['End'] = adata.var['gene_ids'].map(end_mapping)

    return adata


# Take in both 10X feature.tsv and peak_annotation.tsv. Transfer Chromosome, Start, End #
def assign_loc(feature_path, annotation_path, adata):
    
    # Load the full features.tsv and atac_peak_annotation.tsv file 
    features_df = pd.read_csv(feature_path, header=None, sep='\t')
    
    # Set features's header 
    features_df.columns = ['gene_ids', 'gene_symbols', 'feature_type', 'Chromosome', 'Start', 'End']
    features_df.set_index("gene_ids")
    
    # annotation file should already have header 
    annota_df = pd.read_csv(annotation_path, sep='\t')
        

    # merge 2 df by combinations of chrom, start, end
    merged_df = pd.merge(features_df,
                        annota_df,
                        how='left',
                        left_on=['Chromosome', 'Start', 'End'],
                        right_on=['chrom', 'start', 'end']
                        )

    # apply to only atac data
    merged_df = merged_df[merged_df["feature_type"] == "Peaks"]
    
    # make unique by filter rows where 'distance' is 0 and 'gene_id' is a duplicate 
    unique_combinations = merged_df['gene_ids'][merged_df['gene_ids'].duplicated(keep=False)]
    filtered_df = merged_df[(merged_df['gene_ids'].isin(unique_combinations)) & (merged_df['distance'] == 0)]
    
    # when both distance is 0, drop the first one gene 
    unique_combinations = merged_df['gene_ids'].drop_duplicates(keep='last')
    filtered_df = merged_df.loc[unique_combinations.index]

    # Set index as gene_ids
    filtered_df.set_index('gene_ids', inplace=True)

    adata.var['gene'] =  adata.var['gene_ids'].map(filtered_df['gene']) 
    adata.var['distance'] = adata.var['gene_ids'].map(filtered_df['distance']) 
    adata.var['peak_type'] = adata.var['gene_ids'].map(filtered_df['peak_type']) 



    return adata


# separate gene, promotor, CRE for atac data
def separate_GRE_gene_promotor(atac,asisgn_peak_name = 'peak_category', peak = "peak_type",distance = 'distance'):
    
    
    # first separate gene + promotor into gene, other peaks into CRE
    atac.var[asisgn_peak_name] = pd.Series([""] * atac.var.shape[0])
    atac.var[asisgn_peak_name] = atac.var[distance].apply(lambda x: 'gene' if x == 0.0 else 'CRE').copy()
    
    # Create a boolean mask based on the condition
    cre_mask = atac.var[asisgn_peak_name] == 'CRE'
    gene_mask = atac.var[asisgn_peak_name] == 'gene'

    # Use the mask to subset the AnnData object
    adata_CRE = atac[:, cre_mask].copy()
    adata_gene = atac[:, gene_mask].copy()
    
    
    # from gene+promotor, 
    promotor_mask = adata_gene.var[peak] == "promoter"
    gene_mask = adata_gene.var[peak] == "distal"

    adata_promotor = adata_gene[:, promotor_mask].copy()
    adata_promotor.var[asisgn_peak_name] = adata_promotor.var[peak].copy()
    adata_gene = adata_gene[:, gene_mask].copy()
    
    return  adata_CRE, adata_gene, adata_promotor


# separate gene, promotor, CRE for atac data
def separate_GRE_gene(atac,asisgn_peak_name = 'peak_category', peak = "peak_type",distance = 'distance'):
    
    
    # first separate gene + promotor into gene, other peaks into CRE
    atac.var[asisgn_peak_name] = pd.Series([""] * atac.var.shape[0])
    atac.var[asisgn_peak_name] = atac.var[distance].apply(lambda x: 'gene' if x == 0.0 else 'CRE')
    
    # Create a boolean mask based on the condition
    cre_mask = atac.var[asisgn_peak_name] == 'CRE'
    gene_mask = atac.var[asisgn_peak_name] == 'gene'

    # Use the mask to subset the AnnData object
    adata_CRE = atac[:, cre_mask]
    adata_gene = atac[:, gene_mask]
    
    return  adata_CRE, adata_gene


# find gene that is both accessible and expressing in each cell #
def define_open_express_gene(adata_rna, adata_atac, rna_key = "var_names", atac_key = "gene"):
    
    adata_CRE, adata_gene, adata_promoter = separate_GRE_gene_promotor(adata_atac)
       
    # Extract gene names from 'adata'
    rna_gene_names = list(adata_rna.var[rna_key])  
    atac_gene_names = list(adata_gene.var[atac_key])  
    promotor_gene_names = list(adata_promoter.var[atac_key])
        
    # Find the intersection of gene names
    overlap_genes = set(atac_gene_names).intersection(set(rna_gene_names),set(promotor_gene_names) )

    # Subset 'adata' to keep only the overlapping genes
    rna_list = adata_rna.var_names[adata_rna.var[rna_key].isin(overlap_genes)]  
    atac_list = adata_gene.var_names[adata_gene.var[atac_key].isin(overlap_genes)] 
    promoter_list =  adata_promoter.var_names[adata_promoter.var[atac_key].isin(overlap_genes)]
    
    # select shared genes
    adata_rna_flitered = adata_rna[:,rna_list].copy()
    adata_atac_gene_filtered = adata_gene[:,atac_list].copy()
    adata_atac_promoter_filtered = adata_promoter[:,promoter_list].copy()
    
    # add one more label to standardize naming
    adata_rna_flitered.var['gene'] = adata_rna_flitered.var_names

        
    return adata_rna_flitered, adata_atac_gene_filtered, adata_CRE, adata_atac_promoter_filtered


# take atac gene body and atac gene promotor, find gene that is activated in each cell
def find_activated_gene(adata_gene, adata_promotor, key_gene = "gene", key_promotor = "gene"):
    
    # first check they have same gene set
    if (set(adata_gene.var[key_gene]) != set(adata_promotor.var[key_promotor])):
        print("do not have the same gene")
        return 
    
    else:
        # convert the adata value into df
        adata_gene_X= pd.DataFrame(adata_gene.X.toarray(), index=adata_gene.obs_names, columns=adata_gene.var[key_gene])
        adata_gene_X = adata_gene_X.groupby(axis=1, level=0).sum()  # sum the repeated gene 
        adata_gene_X= adata_gene_X[sorted(adata_gene_X.columns)]    # sort the gene in same order
        
        
        adata_promotor_X= pd.DataFrame(adata_promotor.X.toarray(), index=adata_promotor.obs_names, columns=adata_promotor.var[key_promotor])
        adata_promotor_X = adata_promotor_X.groupby(axis=1, level=0).sum()     # sum the repeated gene promotor
        adata_promotor_X= adata_promotor_X[sorted(adata_promotor_X.columns)]   # sort the gene in same order

        
        # for each cell find gene that is both open in its gene body and promotor
        df_intersect = adata_promotor_X * adata_gene_X
        df_binarized = (df_intersect > 0).astype(int)
        
        # reorganize the data into adata formate
        X_sparse = sparse.csr_matrix(df_binarized.values)
        adata_atac = ad.AnnData(X=X_sparse, obs=pd.DataFrame(index=df_binarized.index), var=pd.DataFrame(index=df_binarized.columns))
        
        var_df = df_gene_chromosome(adata_gene) # gather orginal var information
        
        adata_atac.var["Chromosome"] = var_df["Chromosome"] # pass back
        adata_atac.var["Start"] = var_df["Start"]
        adata_atac.var["End"] = var_df["End"]
        adata_atac.var['gene'] = var_df["gene"]
        adata_atac.var['gene_ids'] = var_df["gene_ids"]
        
        
        return adata_atac


# ---------------------------------------------ATAC-preprocessing Functions------------------------------------------------


# ---------------------------------------------RNA-preprocessing Functions------------------------------------------------

# given rna adata, compute average for each cell with  "n_counts" and "n_genes"
# for each gene in each cell, if gene mRNA count > average, save orginal values
# if mRNA count < avereage, save 0
def define_above_baseline_gene(adata):
    
    # Calculate n_counts: total counts per cell (sum of all gene counts per cell)
    adata.obs['n_counts'] = adata.X.sum(axis=1)

    # Calculate n_genes: number of non-zero gene counts per cell (number of expressed genes)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    
    # computer average
    average_mrna_counts = adata.obs['n_counts']/adata.obs['n_genes']
    avg_mrna_per_cell = average_mrna_counts.values

    # perform gene baseline limitation
    expression_matrix = adata.X.toarray()
    new_expression_matrix = np.where(expression_matrix > avg_mrna_per_cell[:, None], expression_matrix, 0)
    
    new_adata = ad.AnnData(X=new_expression_matrix, obs=adata.obs.copy(), var=adata.var.copy())
    
    return new_adata



'''
# limit adata to the shared genes 
def define_open_express_gene(adata_rna, adata_atac, rna_key, atac_key):
    
    adata_CRE, adata_gene = separate_GRE_gene(adata_atac)
       
    # Extract gene names from 'adata'
    rna_gene_names = list(adata_rna.var[rna_key])  
    atac_gene_names = list(adata_gene.var[atac_key])  
        
    # Find the intersection of gene names
    overlap_genes = set(atac_gene_names).intersection(set(rna_gene_names))

    # Subset 'adata' to keep only the overlapping genes
    rna_list = adata_rna.var_names[adata_rna.var[rna_key].isin(overlap_genes)]  
    atac_list = adata_gene.var_names[adata_gene.var[atac_key].isin(overlap_genes)] 
    
    # select genes
    adata_rna_flitered = adata_rna[:,rna_list].copy()
    adata_atac_gene_filtered = adata_gene[:,atac_list].copy()
        
    return adata_rna_flitered, adata_atac_gene_filtered, adata_CRE
'''



# ---------------------------------------------RNA-preprocessing Functions------------------------------------------------


