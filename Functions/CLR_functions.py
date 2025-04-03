### CLR functions
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
import anndata as ad
import pySingleCellNet as pySCN # pip install git+https://github.com/pcahan1/PySingleCellNet.
import random
random.seed(10)
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import normalized_mutual_info_score
from sklearn.feature_selection import mutual_info_regression
from sklearn.feature_selection import mutual_info_classif
from scipy import stats


# isolate gene by chromosome. 
# -1 is for chrX, -2 is for chrY, any positive number is for chr#. # 
def subset_adata_by_chromosome(adata, chromosome_num):
    
    # Convert the chromosome number to its corresponding chromosome label
    if chromosome_num == -1:
        chromosome_label = 'chrX'
    elif chromosome_num == -2:
        chromosome_label = 'chrY'
    else:
        chromosome_label = f'chr{chromosome_num}'

    # Find all features (genes) that are on the specified chromosome
    genes_on_chromosome = adata.var[adata.var['Chromosome'] == chromosome_label].index

    # Subset the adata to include only the columns (genes) that are on the specified chromosome
    adata_subset = adata[:, genes_on_chromosome].copy()

    return adata_subset


# normalize one chromosome's data
def normalize_by_chromasome(adata_rna, original_adata):
    
    # extract selected gene and cells 
    original_adata = original_adata[original_adata.obs.index.isin(adata_rna.obs.index)]
    original_adata = original_adata[:,original_adata.var.index.isin(adata_rna.var.index)]

    # transfer var 
    new_adata = ad.AnnData(X=original_adata.X.copy(), obs=adata_rna.obs.copy(), var=adata_rna.var.copy())

    # Normalization
    sc.pp.normalize_total(new_adata, target_sum=1e4)
    sc.pp.log1p(new_adata)
    
    return new_adata


# normalize all chromsome's data
def normalize_all_chromasome(adata_rna, orginal_count_path, chr_num, chr_y = False):
    
    # get orginal mRNA counts 
    original_adata = sc.read_10x_mtx(orginal_count_path, gex_only = False)
        
    # separate RNA and ATAC data 
    gex_rows = list(map(lambda x: x == 'Gene Expression', original_adata.var['feature_types']))
    original_adata = original_adata[:, gex_rows]
    
    # make variables unique
    original_adata.var_names_make_unique()
    
    # make a new dict to store result
    normalized_all_chr = {}
    
    # normalize by chr # 
    for i in range(1, chr_num + 1):
        chromosome_subsets = subset_adata_by_chromosome(adata_rna,i)
        new_adata = normalize_by_chromasome(chromosome_subsets,original_adata)
        
        normalized_all_chr[f'chr{i}'] = new_adata
        
    # normalize chrX
    chromosome_subsets = subset_adata_by_chromosome(adata_rna,-1)
    new_adata = normalize_by_chromasome(chromosome_subsets,original_adata)
    normalized_all_chr['chrX'] = new_adata
    
    # normalize chrY
    if (chr_y == True):
        chromosome_subsets = subset_adata_by_chromosome(adata_rna,-2)
        new_adata = normalize_by_chromasome(chromosome_subsets,original_adata)
        normalized_all_chr['chrY'] = new_adata
    
    return normalized_all_chr


# Calculate the Mutual Information score for gene expressions and chromatin accessbility 
def MI_Matrix(adata_atac, adata_rna):
    
    if adata_atac.n_obs != adata_rna.n_obs:
        print("The two datasets do not have the same number of cells.")
        return  # Exit the function early
    
    # convert rna and atac adata object to dataframe 
    rna_expression = adata_rna.X.toarray()
    rna_df = pd.DataFrame(rna_expression, index=adata_rna.obs_names, columns=adata_rna.var_names)

    atac_access = adata_atac.X.toarray()
    atac_df = pd.DataFrame(atac_access, index=adata_atac.obs_names, columns=adata_atac.var_names)
    
    # normalize by KBinsDiscretizer
    est=KBinsDiscretizer(n_bins=int(np.sqrt(atac_df.shape[1])),encode='ordinal',strategy="quantile")
    est.fit(atac_df)
    dataset_atac=est.transform(atac_df)

    est=KBinsDiscretizer(n_bins=int(np.sqrt(rna_df.shape[1])),encode='ordinal',strategy="quantile")
    est.fit(rna_df)
    dataset_rna=est.transform(rna_df)
    
    # calculate MI matrix 
    n_cols_rna = dataset_rna.shape[1]  #number of gene
    n_cols_atac = dataset_atac.shape[1]  #number of peak

    # Initialize the mutual information matrix
    mi_matrix = np.zeros((n_cols_rna, n_cols_atac))

    # Calculate mutual information for each pair of columns
    for i in range(n_cols_rna):
        for j in range(n_cols_atac):
                mi_matrix[i, j] = mutual_info_score(dataset_rna[:, i], dataset_atac[:, j])
               
    # convert to df  
    gene_names = adata_rna.var_names
    peak_names = adata_atac.var["gene_ids"]
                
    return mi_matrix


# Calculate the Mutual Information score for gene expressions and chromatin accessbility 
def MI_Matrix_MIinfoClassif(adata_atac, adata_rna):
    
    if adata_atac.n_obs != adata_rna.n_obs:
        print("The two datasets do not have the same number of cells.")
        return  # Exit the function early
    
    # convert rna and atac adata object to dataframe 
    rna_expression = adata_rna.X.toarray()
    atac_access = adata_atac.X.toarray() 
    atac_access = (atac_access > 0).astype(int)   
    
    # calculate MI matrix 
    n_cols_rna = rna_expression.shape[1]  #number of gene
    n_cols_atac = atac_access.shape[1]  #number of peak

    # Initialize the mutual information matrix
    mi_matrix = np.zeros((n_cols_rna, n_cols_atac))

    # Calculate mutual information for each pair of columns
    for i in range(n_cols_rna):
        rna = rna_expression[:,i].reshape(-1, 1)
        for j in range(n_cols_atac):
                atac = atac_access[:, j]
                mi = mutual_info_classif(rna, atac, discrete_features = False)
                mi_matrix[i, j] = mi[0]
               
    # convert to df  
    mi_df = pd.DataFrame(mi_matrix, index=adata_rna.var_names, columns=adata_atac.var_names)
                
    return mi_df


# Computer z-score matrix from MI matrix
def CLR_Matrix(mi_matrix):
    
    # calculate z-score
    z_score_atac=pd.DataFrame(stats.zscore(mi_matrix,axis=0))
    z_score_atac[z_score_atac<0]=0

    z_score_rna=pd.DataFrame(stats.zscore(mi_matrix,axis=1))
    z_score_rna[z_score_rna<0]=0

    # combine z-score
    CLR_Matrix =np.sqrt(z_score_atac**2+z_score_rna**2)

    return CLR_Matrix

"""
# start with ChrY, ChrX, then Chr1 to the rest
def chromosome_correlation(chromosom_num, adata_rna, adata_atac, chromosome_X = True, chromosome_Y = True):
    
    MI = {}
    
    for i in range (1,chromosom_num):
        rna_chrom_sub = subset_adata_by_chromosome(adata_rna, i)
        atac_chrom_sub = subset_adata_by_chromosome(adata_atac, i)
        
        chromosome_subsets = MI_Matrix(rna_chrom_sub,atac_chrom_sub)
        MI[f'chr{chromosom_num}'] = chromosome_subsets
    

    return MI
    
"""


"""
# Compute CLR matrix for Gene peaks with non-gene open peaks 
def MI_CREtoGene_by_Chromosome(adata_rna, adata_atac,chr_num, gene_list = None, rna_key = None, atac_key = None):
    
    adata_rna, adata_atac_gene, adata_atac_CRE = define_open_express_gene(adata_rna,adata_atac,"var_names", "gene")
    
    
    chr_gem = subset_adata_by_chromosome(adata_rna,chr_num)
    chr_atac_gene = subset_adata_by_chromosome(adata_atac_gene,chr_num) 
    chr_atac_cre = subset_adata_by_chromosome(adata_atac_CRE,chr_num) 
    
    if (gene_list is not None) & (atac_key is not None) & (rna_key is not None) :
          chr_gem = limit_gene_by_list(gene_list, chr_gem, rna_key)
          chr_atac_gene = limit_gene_by_list(gene_list, chr_atac_gene, atac_key)
    
    
    MI_peaks = MI_Matrix(chr_atac_gene,chr_atac_cre)
    CLR_peaks = CLR_Matrix(MI_peaks)
    
    MI_rna = MI_Matrix(chr_gem,chr_atac_cre)
    CLR_rna = CLR_Matrix(MI_rna)
    
    if (rna_key is not None)& (atac_key is not None):
          CLR_peaks.index = list(chr_atac_cre.var_names)
          CLR_peaks.columns =  list(chr_atac_gene.var[atac_key])
          
          CLR_rna.index =  list(chr_atac_cre.var_names)
          CLR_rna.columns = list(chr_gem.var[rna_key]) 
     

    
    return CLR_peaks, CLR_rna, MI_peaks, MI_rna

"""


"""
# Compute CLR matrix for Gene expresses with non-gene open peaks
def MI_TFtoCRE_by_Chromosome(adata_rna, adata_atac,chr_num, TF_file = None, rna_key = None, atac_key = None):
    
    adata_rna, adata_atac_gene, adata_atac_CRE = define_open_express_gene(adata_rna,adata_atac,"var_names", "gene")
    
    chr_gem = subset_adata_by_chromosome(adata_rna,chr_num)
    chr_atac_cre = subset_adata_by_chromosome(adata_atac_CRE,chr_num) 
    
    if (TF_file is not None) & (rna_key is not None) :
          chr_gem = limit_gene_by_file(TF_file, chr_gem, rna_key)
          print(chr_gem)
        
    #MI_rna = MI_Matrix(chr_gem,chr_atac_cre)
    #CLR_rna = CLR_Matrix(MI_rna)
    
    #if (rna_key is not None)& (atac_key is not None):
     #     CLR_rna.index =  list(chr_atac_cre.var_names)
      #    CLR_rna.columns = list(chr_gem.var[rna_key]) 
          
    
    return chr_gem#CLR_rna

"""


"""
def filter_by_zscore(CLR_peaks,CLR_rna, zscore):
    CLR_peaks[CLR_peaks < zscore] = 0
    CLR_rna[CLR_rna < zscore] = 0
    return CLR_peaks.copy(), CLR_rna.copy()
"""



# takes the two adata set with same cells to perform MI association between CRE peaks and gene peaks
# gene_names and cre_names are names for CRE peaks and gene stored in orginal data
def compute_CLR_CRE_chromosome(adata_cre, adata_gene, start_chr , end_chr ,  gene_names = 'gene', cre_names = 'gene_ids', chr_x = False, chr_y = False):
    
    MI_matrixes = {}
    
    # run CLR on X chr
    def Chr_x():
        
        # subset 
        adata_gene_sub = subset_adata_by_chromosome(adata_gene, -1)
        adata_cre_sub = subset_adata_by_chromosome(adata_cre, -1)
        
        # compute MI 
        chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
        
        # rename 
        chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])
        
        MI_matrixes[f'chrX'] = chromosome_subsets_df
        print('completed chrX')
        
    # run CLR on y chr
    def Chr_y():
        adata_gene_sub = subset_adata_by_chromosome(adata_gene, -2)
        adata_cre_sub = subset_adata_by_chromosome(adata_cre, -2)
            
        chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
        chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])
        
        MI_matrixes['chrY'] = chromosome_subsets_df
        print('completed chrY')
        
    # run CLR on # chr
    def chr_n():
        for i in range (start_chr,end_chr+1):
            adata_gene_sub = subset_adata_by_chromosome(adata_gene, i)
            adata_cre_sub = subset_adata_by_chromosome(adata_cre, i)
            
            chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
            chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])

            MI_matrixes[f'chr{i}'] = chromosome_subsets_df
            print(f'completed chr{i}')
    
    
    # only wants x or y chr
    if ((start_chr <= 0) or (end_chr <= 0)):
        
        if chr_x == True: Chr_x()
        if chr_y == True: Chr_y()
        return MI_matrixes  
    
    # wants # chr with x or y chr 
    else:
        chr_n()
        if chr_x == True: Chr_x()
        if chr_y == True: Chr_y() 
        return MI_matrixes
    

# compute for CLR matrixes  for all chromosomes
def compute_CLR_matrixes(MI_matrixes):
    CLR_matrixes =  {}
    
    # iterate through every chr
    for key, value in MI_matrixes.items():
        result = CLR_Matrix(value) # compute value for each chr 
        CLR_matrixes[key] = result
        
    return CLR_matrixes


# takes the two adata set with same cells to perform MI association between CRE peaks and gene expressions
# gene_names and cre_names are names for CRE peaks and gene stored in orginal data
def compute_CLR_RNA_chromosome(adata_cre, adata_gene, start_chr, end_chr,  gene_names = 'var_names', cre_names = 'gene_ids', chr_x = False, chr_y = False):
    
    
    MI_matrixes = {}
    
    # run CLR on X chr
    def Chr_x():
        
        # subset 
        adata_gene_sub =  adata_gene['chrX']
        adata_cre_sub = subset_adata_by_chromosome(adata_cre, -1)
        
        # compute MI 
        chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
        
        # rename 
        chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])
        
        MI_matrixes['chrX'] = chromosome_subsets_df
        print('completed chrX')
        
    # run CLR on y chr
    def Chr_y():
        adata_gene_sub = adata_gene['chrY'] 
        adata_cre_sub = subset_adata_by_chromosome(adata_cre, -2)
            
        chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
        chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])
        
        MI_matrixes['chrY'] = chromosome_subsets_df
        print('completed chrY')
        
    # run CLR on # chr
    def chr_n():
        for i in range (start_chr,end_chr+1):
            adata_gene_sub =  adata_gene[f'chr{i}']
            adata_cre_sub = subset_adata_by_chromosome(adata_cre, i)
            
            chromosome_subsets = MI_Matrix(adata_cre_sub,adata_gene_sub)
            chromosome_subsets_df = pd.DataFrame(chromosome_subsets, columns=adata_cre_sub.var[cre_names], index = adata_gene_sub.var[gene_names])

            MI_matrixes[f'chr{i}'] = chromosome_subsets_df
            print(f'completed chr{i}')
    
    
    # only wants x or y chr
    if (start_chr <= 0 or end_chr <= 0):
        
        if chr_x == True: Chr_x()
        if chr_y == True: Chr_y()
        return MI_matrixes  
    
    # wants # chr with x or y chr 
    else:
        chr_n()
        if chr_x == True: Chr_x()
        if chr_y == True: Chr_y() 
        return MI_matrixes





# filter all CLR_matrixes by z-score
def z_score_filter(CLR_matrixes, zthre):
    zthre_matriex = {}
    
        # iterate through every chr
    for key, value in CLR_matrixes.items():
        
        # filter out all entries smaller than z-thre
        result = value.applymap(lambda x: x if x > zthre else 0 )
        
        # drop any gene or peaks that has no association to ther peaks or genes
        result = result.loc[:,value.sum(axis=0) != 0]
        result = result.loc[value.sum(axis=1) != 0]
        
        zthre_matriex[key] = result
        
    return zthre_matriex


# intesect two CLR matrixes 
def intersect_chromsome_matrixes(CLR_1, CLR_2, chr_x = False, chr_y = False):
    
    intersected_chr_matrixes = {}
    
    # for chr #
    for key in CLR_1:
        if key in CLR_2:
            chr_CLR_1 = CLR_1[key]
            chr_CLR_2 = CLR_2[key]
            
            chr_CLR_1 = chr_CLR_1.sort_index()
            chr_CLR_2 = chr_CLR_2.sort_index()
            
            intersected_chr_matrixes[key] = chr_CLR_1*chr_CLR_2
            
    # for chrX
    if chr_x == True:
        chr_CLR_1 = CLR_1['chrX']
        chr_CLR_2 = CLR_2['chrX']
        
        chr_CLR_1 = chr_CLR_1.sort_index()
        chr_CLR_2 = chr_CLR_2.sort_index()
        
        intersected_chr_matrixes['chrX'] = chr_CLR_1*chr_CLR_2

    # for chrY
    if chr_y == True:
        chr_CLR_1 = CLR_1['chrY']
        chr_CLR_2 = CLR_2['chrY']
        
        chr_CLR_1 = chr_CLR_1.sort_index()
        chr_CLR_2 = chr_CLR_2.sort_index()
        
        intersected_chr_matrixes['chrY'] = chr_CLR_1*chr_CLR_2
    
    return intersected_chr_matrixes





"""

# intersect 
def CLR_intersect(CLR_peaks,CLR_rna):
    
    CLR_peaks = CLR_peaks.fillna(0.00)
    CLR_rna = CLR_rna.fillna(0.00)

    CLR_peaks = CLR_peaks.groupby(axis=1, level=0).sum()
    
    CLR_peaks = CLR_peaks.reindex(sorted(CLR_peaks.columns), axis=1)
    CLR_rna = CLR_rna.reindex(sorted(CLR_rna.columns), axis=1)

    CLR = CLR_peaks + CLR_rna 
    
    non_cor_peaks = CLR.sum(axis=1)

    # Drop rows where the sum is zero
    empty_peaks = CLR[non_cor_peaks == 0]  
    peaks = CLR[non_cor_peaks != 0]  


    return CLR, empty_peaks, peaks

    Returns:
        _type_: _description_
    """

'''
# intersect 
def CLR_intersect(CLR_peaks,CLR_rna):
    
    CLR_peaks = CLR_peaks.fillna(0.00)
    CLR_rna = CLR_rna.fillna(0.00)

    CLR_peaks = CLR_peaks.groupby(axis=1, level=0).sum()
    
    CLR_peaks = CLR_peaks.reindex(sorted(CLR_peaks.columns), axis=1)
    CLR_rna = CLR_rna.reindex(sorted(CLR_rna.columns), axis=1)

    CLR = CLR_peaks + CLR_rna 
    
    non_cor_peaks = CLR.sum(axis=1)

    # Drop rows where the sum is zero
    empty_peaks = CLR[non_cor_peaks == 0]  
    peaks = CLR[non_cor_peaks != 0]  


    return CLR, empty_peaks, peaks

# start with ChrY, ChrX, then Chr1 to the rest
def chromosome_correlation(chromosom_num, adata_rna, adata_atac, chromosome_X = True, chromosome_Y = True):
    
    MI = {}
    
    for i in range (1,chromosom_num):
        rna_chrom_sub = subset_adata_by_chromosome(adata_rna, i)
        atac_chrom_sub = subset_adata_by_chromosome(adata_atac, i)
        
        chromosome_subsets = MI_Matrix(rna_chrom_sub,atac_chrom_sub)
        MI[f'chr{chromosom_num}'] = chromosome_subsets
    

    return MI

# Compute CLR matrix for Gene peaks with non-gene open peaks 
def MI_CREtoGene_by_Chromosome(adata_rna, adata_atac,chr_num, gene_list = None, rna_key = None, atac_key = None):
    
    adata_rna, adata_atac_gene, adata_atac_CRE = define_open_express_gene(adata_rna,adata_atac,"var_names", "gene")
    
    
    chr_gem = subset_adata_by_chromosome(adata_rna,chr_num)
    chr_atac_gene = subset_adata_by_chromosome(adata_atac_gene,chr_num) 
    chr_atac_cre = subset_adata_by_chromosome(adata_atac_CRE,chr_num) 
    
    if (gene_list is not None) & (atac_key is not None) & (rna_key is not None) :
          chr_gem = limit_gene_by_list(gene_list, chr_gem, rna_key)
          chr_atac_gene = limit_gene_by_list(gene_list, chr_atac_gene, atac_key)
    
    
    MI_peaks = MI_Matrix(chr_atac_gene,chr_atac_cre)
    CLR_peaks = CLR_Matrix(MI_peaks)
    
    MI_rna = MI_Matrix(chr_gem,chr_atac_cre)
    CLR_rna = CLR_Matrix(MI_rna)
    
    if (rna_key is not None)& (atac_key is not None):
          CLR_peaks.index = list(chr_atac_cre.var_names)
          CLR_peaks.columns =  list(chr_atac_gene.var[atac_key])
          
          CLR_rna.index =  list(chr_atac_cre.var_names)
          CLR_rna.columns = list(chr_gem.var[rna_key]) 
     

    
    return CLR_peaks, CLR_rna, MI_peaks, MI_rna

# Compute CLR matrix for Gene expresses with non-gene open peaks
def MI_TFtoCRE_by_Chromosome(adata_rna, adata_atac,chr_num, TF_file = None, rna_key = None, atac_key = None):
    
    adata_rna, adata_atac_gene, adata_atac_CRE = define_open_express_gene(adata_rna,adata_atac,"var_names", "gene")
    
    chr_gem = subset_adata_by_chromosome(adata_rna,chr_num)
    chr_atac_cre = subset_adata_by_chromosome(adata_atac_CRE,chr_num) 
    
    if (TF_file is not None) & (rna_key is not None) :
          chr_gem = limit_gene_by_file(TF_file, chr_gem, rna_key)
          print(chr_gem)
        
    #MI_rna = MI_Matrix(chr_gem,chr_atac_cre)
    #CLR_rna = CLR_Matrix(MI_rna)
    
    #if (rna_key is not None)& (atac_key is not None):
     #     CLR_rna.index =  list(chr_atac_cre.var_names)
      #    CLR_rna.columns = list(chr_gem.var[rna_key]) 
          
    
    return chr_gem#CLR_rna

'''
