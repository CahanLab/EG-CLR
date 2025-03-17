import pySingleCellNet as pySCN
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import normalized_mutual_info_score
from scipy import stats

from .QC_functions import *



### CLR functions

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


# isolate gene by chromosome. 
# -1 is for chrX, -2 is for chrY, any positive number is for chr#. 
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
    est=KBinsDiscretizer(n_bins=int(np.sqrt(atac_df.shape[1])),encode='ordinal',strategy="uniform")
    est.fit(atac_df)
    dataset_atac=est.transform(atac_df)

    est=KBinsDiscretizer(n_bins=int(np.sqrt(rna_df.shape[1])),encode='ordinal',strategy="uniform")
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
                mi_matrix[i, j] = normalized_mutual_info_score(dataset_rna[:, i], dataset_atac[:, j])
               
    # convert to df  
    gene_names = adata_rna.var_names
    peak_names = adata_atac.var["gene_ids"]
                
    return mi_matrix

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


def filter_by_zscore(CLR_peaks,CLR_rna, zscore):
    CLR_peaks[CLR_peaks < zscore] = 0
    CLR_rna[CLR_rna < zscore] = 0
    return CLR_peaks.copy(), CLR_rna.copy()

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

