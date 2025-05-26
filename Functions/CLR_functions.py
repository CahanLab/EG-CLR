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
from sklearn.metrics import mutual_info_score
from sklearn.feature_selection import mutual_info_regression
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import StandardScaler
from scipy import stats
from joblib import Parallel, delayed



# ---------------------------------------------CLR Data Preprocessing Functions------------------------------------------------


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


# replace one chromosome's data by original counts
def replace_by_origin(adata, original_adata):
    
    # extract selected gene and cells 
    original_adata = original_adata[original_adata.obs.index.isin(adata.obs.index)]
    original_adata = original_adata[:,original_adata.var.index.isin(adata.var.index)]

    # transfer var 
    new_adata = ad.AnnData(X=original_adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())

    
    return new_adata


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



# ---------------------------------------------CLR Data Preprocessing Functions------------------------------------------------






# ---------------------------------------------Gene filtering Functions------------------------------------------------


# separate gene, promotor, CRE for atac data #
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


# separate gene, promotor + CRE for atac data #
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


# find gene that is both accessible in gene body and promoter, and expressing in each cell #
def define_rna_promoter_gene_OCregion(adata_rna, adata_atac, atac_key = "gene"):
    
    adata_CRE, adata_gene, adata_promoter = separate_GRE_gene_promotor(adata_atac)
           
    # Extract gene names from 'adata'
    rna_gene_names = list(adata_rna.var_names)  
    atac_gene_names = list(adata_gene.var[atac_key])  
    promotor_gene_names = list(adata_promoter.var[atac_key])
        
    # Find the intersection of gene names
    overlap_genes = set(atac_gene_names).intersection(set(rna_gene_names),set(promotor_gene_names) )

    # Subset 'adata' to keep only the overlapping genes
    rna_list = adata_rna.var_names[adata_rna.var_names.isin(overlap_genes)]  
    atac_list = adata_gene.var_names[adata_gene.var[atac_key].isin(overlap_genes)] 
    promoter_list =  adata_promoter.var_names[adata_promoter.var[atac_key].isin(overlap_genes)]
    
    # select shared genes
    adata_rna_flitered = adata_rna[:,rna_list].copy()
    adata_atac_gene_filtered = adata_gene[:,atac_list].copy()
    adata_atac_promoter_filtered = adata_promoter[:,promoter_list].copy()
    
    # add one more label to standardize naming
    adata_rna_flitered.var['gene'] = adata_rna_flitered.var_names

        
    return adata_rna_flitered, adata_atac_gene_filtered, adata_CRE, adata_atac_promoter_filtered


# find gene that is both accessible in promoter, and expressing in each cell #
def define_rna_promoter_OCregion(adata_rna, adata_atac, atac_key = "gene"):
    
    _, _, adata_promoter = separate_GRE_gene_promotor(adata_atac)
    
    # first check they have same gene set
    if (set(adata_rna.var_names).intersection(set(adata_promoter.var["gene"])) == 0):
        print("do not have the same gene")
        return 
           
    # Extract gene names from 'adata'
    rna_gene_names = list(adata_rna.var_names)  
    promotor_gene_names = list(adata_promoter.var[atac_key])
        
    # Find the intersection of gene names
    overlap_genes = set(promotor_gene_names).intersection(set(rna_gene_names))

    # Subset 'adata' to keep only the overlapping genes
    rna_list = adata_rna.var_names[adata_rna.var_names.isin(overlap_genes)]  
    promoter_list =  adata_promoter.var_names[adata_promoter.var[atac_key].isin(overlap_genes)]
    
    # select shared genes
    adata_rna_flitered = adata_rna[:,rna_list].copy()
    adata_atac_promoter_filtered = adata_promoter[:,promoter_list].copy()
    
    # add one more label to standardize naming
    adata_rna_flitered.var['gene'] = adata_rna_flitered.var_names

        
    return adata_rna_flitered, adata_atac_promoter_filtered


# combine promoter and gene body counts #
def combine_promoter_gene_counts(adata_promoter, adata_gene_body, gene_key="gene"):
    
    adata = ad.concat(
        [adata_promoter, adata_gene_body],
        axis         = 1,           
        join         = "outer",    
        merge        = "same",      
        label        = None,        
        index_unique = None         
    )
    
    adata_df = adata.to_df().transpose()
    adata_df[gene_key] = adata.var[gene_key]
    
    gene_sums = adata_df.groupby(gene_key).sum(numeric_only=True)
    gene_sums_sparse= sparse.csr_matrix(gene_sums.transpose().values)  
    
    return ad.AnnData(X=gene_sums_sparse, obs=pd.DataFrame(index=gene_sums.columns), var = pd.DataFrame(index=gene_sums.index))


# ---------------------------------------------Gene filtering Functions------------------------------------------------





# ---------------------------------------------CLR Core Functions---------------------------------------------------------


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


# calculate mutual information between RNA and ATAC data. 
# Feature matrix is continous gene expression data. Target column is binary chromatin accessbility data.
# n_neighbors is the number of neighbors to use for KNN estimator.
# n_jobs is the number of jobs to run in parallel. #
def MI_Matrix_MIinfoClassif( adata_rna, adata_atac, *, n_neighbors=3, n_jobs=1, random_state=0):

    if adata_atac.n_obs != adata_rna.n_obs:
        print("The two datasets do not have the same number of cells.")
        return  # Exit the function early
    
    # convert rna and atac adata object to dataframe 
    rna_expression = adata_rna.X.toarray()
    atac_access = adata_atac.X.toarray() 
    
    # Use KNN entropy estimator which works in Euclidean space. Every variable in feature should be on the 
    # same scale to avoid large numeric ranges dominate the neighbor
    # This can be done by give every gene x-mu/sigma (so variance is 1)
    scaler = StandardScaler()              # mean‑0, std‑1 per column
    rna_expression_scaled = scaler.fit_transform(rna_expression)     # X is (n_samples, n_features)   
    atac_access_scaled = scaler.fit_transform(atac_access)     # X is (n_samples, n_features)   
     
    atac_access_scaled = (atac_access_scaled > 0).astype(int) 

    # helper function to compute mutual information
    def _mi_info(j) -> np.ndarray:
        y = atac_access_scaled[:, j]
        return (mutual_info_classif(rna_expression_scaled, y, n_neighbors=n_neighbors, n_jobs=n_jobs, random_state=random_state, discrete_features=False))
        
    # parallelize
    mi_list = Parallel(n_jobs=n_jobs)(delayed(_mi_info)(j) for j in range(atac_access_scaled.shape[1]))
    
    mi_matrix = np.vstack(mi_list).T
               
    # convert to df  
    mi_df = pd.DataFrame(mi_matrix, index=adata_atac.var_names, columns=adata_rna.var_names)
                
    return mi_df


# calculate mutual information between RNA and ATAC data. 
# Feature matrix is continous gene expression data. Target column is continous chromatin accessbility data.
# n_neighbors is the number of neighbors to use for KNN estimator.
# n_jobs is the number of jobs to run in parallel. #
def MI_Matrix_MIinfoRegression(adata_rna, adata_atac, *, n_neighbors=3, n_jobs=1, random_state=0):
    
    if adata_atac.n_obs != adata_rna.n_obs:
        print("The two datasets do not have the same number of cells.")
        return  # Exit the function early
    
    # convert rna and atac adata object to dataframe 
    rna_expression = adata_rna.X.toarray()
    atac_access = adata_atac.X.toarray()
    
    
    # Use KNN entropy estimator which works in Euclidean space. Every variable in feature should be on the 
    # same scale to avoid large numeric ranges dominate the neighbor
    # This can be done by give every gene x-mu/sigma (so variance is 1)
    scaler = StandardScaler()              # mean‑0, std‑1 per column
    rna_expression_scaled = scaler.fit_transform(rna_expression)     # X is (n_samples, n_features)   
    atac_access_scaled = scaler.fit_transform(atac_access)     # X is (n_samples, n_features)
     
     
    # helper function to compute mutual information
    def _mi_info(j) -> np.ndarray:
        y = atac_access_scaled[:, j]
        return (mutual_info_regression(rna_expression_scaled, y, n_neighbors=n_neighbors, n_jobs=n_jobs, random_state=random_state, discrete_features=False))
        
    # parallelize
    mi_list = Parallel(n_jobs=n_jobs)(delayed(_mi_info)(j) for j in range(atac_access_scaled.shape[1]))

    
    mi_matrix = np.vstack(mi_list).T
               
    # convert to df  
    mi_df = pd.DataFrame(mi_matrix, index=adata_atac.var_names, columns=adata_rna.var_names)
                
    return mi_df
    
    if adata_atac.n_obs != adata_rna.n_obs:
        print("The two datasets do not have the same number of cells.")
        return  # Exit the function early
    
    # convert rna and atac adata object to dataframe 
    rna_expression = adata_rna.X.toarray()
    atac_access = adata_atac.X.toarray() 
    
    
    '''
    # Use KNN entropy estimator which works in Euclidean space. Every variable in feature should be on the 
    # same scale to avoid large numeric ranges dominate the neighbor
    # This can be done by give every gene x-mu/sigma (so variance is 1)
    scaler = StandardScaler()              # mean‑0, std‑1 per column
    rna_expression_scaled = scaler.fit_transform(rna_expression)     # X is (n_samples, n_features)   
    '''
     
    # helper function to compute mutual information
    def _mi_info(j) -> np.ndarray:
        y = atac_access[:, j]
        return (mutual_info_regression(rna_expression, y, n_neighbors=n_neighbors, n_jobs=n_jobs, random_state=random_state, discrete_features=False))
        
    # parallelize
    mi_list = Parallel(n_jobs=n_jobs)(delayed(_mi_info)(j) for j in range(atac_access.shape[1]))
    
    mi_matrix = np.vstack(mi_list).T
               
    # convert to df  
    mi_df = pd.DataFrame(mi_matrix, index=adata_atac.var_names, columns=adata_rna.var_names)
                
    return mi_df


# Given the MI matrix, adjust the MI values by peak distance to target gene TSS
def distance_adjustment(MI_Matrix, gene_tss):
    
    MI_Matrix_adjusted = MI_Matrix.copy()

    for gene_name in MI_Matrix.columns:
        
        gene = MI_Matrix[gene_name] # obtain peaks name and MI values assciated with the gene
        tss = gene_tss[gene_name] # obtain TSS of the gene

        for peak_name,value in gene.items():
            
            # for each peak, obtain the distance from middle of peak to TSS
            start = peak_name.split(':')[1].split('-')[0]
            end = peak_name.split(':')[1].split('-')[1]
            mid = int((int(start) + int(end)) / 2)
            
            distance = abs(mid - tss)
            
            #print(f"Gene: {gene_name}, tss: {tss}, Peak: {peak_name}, mid : {mid}, Distance to TSS: {distance}")
            
            # Ajust the MI value by distanc
            MI_adjusted = (gene[peak_name]* (1/distance))
            MI_Matrix_adjusted[gene_name][peak_name] = MI_adjusted
        
        
    return MI_Matrix_adjusted
    
    
# Computer z-score matrix from MI matrix #
def CLR_Matrix(mi_matrix):
    
    # calculate z-score
    z_score_atac=pd.DataFrame(stats.zscore(mi_matrix,axis=0,ddof=0, nan_policy='omit'))
    z_score_atac[z_score_atac<0]=0

    z_score_rna=pd.DataFrame(stats.zscore(mi_matrix,axis=1,ddof=0, nan_policy='omit'))
    z_score_rna[z_score_rna<0]=0

    # combine z-score
    CLR_Matrix =np.sqrt(z_score_atac**2+z_score_rna**2)
    
    # rename the columns and index
    CLR_Matrix = pd.DataFrame(CLR_Matrix.values, index=mi_matrix.index, columns=mi_matrix.columns)

    return CLR_Matrix


# filter all CLR_matrixes by z-score
def z_score_filter(CLR_matrix, zthre):
    
    # filter out all entries smaller than z-thre
    result = CLR_matrix.applymap(lambda x: x if x > zthre else 0.0 )
        
    # drop any gene or peaks that has no association to ther peaks or genes
    result = result.loc[:,CLR_matrix.sum(axis=0) != 0]
    result = result.loc[CLR_matrix.sum(axis=1) != 0]
        
    return result


# intesect two CLR matrixes 
def intersect_chromsome_matrixes(CLR_1, CLR_2):
    
    CLR_1 = CLR_1.sort_index(axis=0).sort_index(axis=1) 
    CLR_2 = CLR_2.sort_index(axis=0).sort_index(axis=1) 
            
    intersected_matrix = np.sqrt(CLR_1**2+CLR_2**2)
    
    intersected_matrix = intersected_matrix.loc[:,intersected_matrix.sum(axis=0) != 0]
    intersected_matrix = intersected_matrix.loc[intersected_matrix.sum(axis=1) != 0]
    
    return pd.DataFrame(intersected_matrix, index=CLR_1.index, columns=CLR_1.columns)


# ---------------------------------------------CLR Core Functions---------------------------------------------------------




# ---------------------------------------------Scale up Functions---------------------------------------------------------


# compute for CLR matrixes  for all chromosomes
def compute_CLR_matrixes(MI_matrixes):
    
    CLR_matrixes =  {}
    
    # iterate through every chr
    for key, value in MI_matrixes.items():
        result = CLR_Matrix(value) # compute value for each chr 
        CLR_matrixes[key] = result
        
    return CLR_matrixes


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


# ---------------------------------------------Scale up Functions---------------------------------------------------------




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
