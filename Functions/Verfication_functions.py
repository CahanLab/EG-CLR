### verfication
import pandas as pd
from muon import atac as ac
from muon import prot as pt
import anndata as ad
import random
import numpy as np
random.seed(10)
import pybedtools


# ---------------------------------------------Speed Test Functions------------------------------------------------


def random_subset_adatas(adata1: ad.AnnData, adata2: ad.AnnData, n_cells: int, random_state: int | None = None, copy: bool = True ):

    # sanity checks  
    if adata1.n_obs != adata2.n_obs or not np.array_equal(
        adata1.obs_names, adata2.obs_names
    ):
        raise ValueError("Both AnnData objects must share the same cells in the same order.")
    if n_cells > adata1.n_obs:
        raise ValueError(f"n_cells={n_cells} is greater than total cells ({adata1.n_obs}).")

    # sampling  
    rng = np.random.default_rng(random_state)
    idx = rng.choice(adata1.n_obs, size=n_cells, replace=False)

    ad1_sub = adata1[idx].copy() if copy else adata1[idx]
    ad2_sub = adata2[idx].copy() if copy else adata2[idx]

    return ad1_sub, ad2_sub

# ---------------------------------------------Speed Test Functions------------------------------------------------




# ---------------------------------------------Data Exploration Functions------------------------------------------------


# Turn CLR matrix in to DF with 3 columns: Element, Value, Gene
def convert_matrix_to_df(matrix):
    
    results = []

    for i in range(matrix.shape[1]):
        # Filter positive values and extract the column as Series
        prediction_list = matrix[matrix.iloc[:, i] > 0.0].iloc[:, i].copy()
        
        # Get the name of the column
        name = matrix.columns[i]
        
        # Append each (index, value, name) triplet as a dict
        for idx, val in prediction_list.items():
            results.append({'Element': idx, 'Value': val, 'Gene': name})

    # Create final DataFrame
    df_result = pd.DataFrame(results)

    return df_result
    

def calculate_tss(prediction_df, adata):
    
    # extract gene coordinate
    adata.var["Gene"] = adata.var_names
    gene_coordinate = adata.var[["Gene", "Chromosome", "Start", "End"]]

    # add gene coordinate to CRISPRi data
    prediction_df_tss = (
        prediction_df.merge(gene_coordinate[['Gene',"Start", "End"]],        
            left_on='Gene',              
            right_on='Gene',            
            how='left',
            validate="many_to_one")                  
     )

    # expand peak coordinate
    prediction_df_tss[['Chr_peak', 'Start_peak', 'End_peak']] = prediction_df_tss['Element'].str.extract(r'(chr[^:]+):(\d+)-(\d+)')
    prediction_df_tss['Start_peak'] = prediction_df_tss['Start_peak'].astype(int)
    prediction_df_tss['End_peak'] = prediction_df_tss['End_peak'].astype(int)

    # calcaulte tss (middle of gene coordiante to middle of peak coordinate)
    prediction_df_tss["distance_to_tss"] = ((prediction_df_tss["Start"] + prediction_df_tss["End"]) // 2  
                                            - (prediction_df_tss["Start_peak"] + prediction_df_tss["End_peak"] // 2)).abs()

    return prediction_df_tss


# ---------------------------------------------Data Exploration  Functions------------------------------------------------




# ---------------------------------------------CRRSPRi Functions------------------------------------------------


# Take the prediction and CRIPSRi data for a gene (both must contain chr, start, end)
# ensure ATAC-seq data and CRISPRi data are in the same genomic range #
def select_testing_peaks(CRISPRi, ATAC):
    
    # Find the range of CRISPRi data
    CRISPRi_start =  CRISPRi['start'].min()
    CRISPRi_end = CRISPRi['end'].max() 

    #print("CRISPRi starts at ",CRISPRi_start, " end at ", CRISPRi_end)

    # Filter for overlapping peaks
    ATAC_limited = ATAC[
        (ATAC["end"].values.astype(int) <  CRISPRi_end) & 
        (ATAC["start"].values.astype(int) > CRISPRi_start)
    ].copy()
    
    ATAC_limited[["start", "end"]] = ATAC_limited[["start", "end"]].astype(int)

    return ATAC_limited


# coverlap between two data sets (chr, start, end)#
def CRISPRi_comparison(A, B):
    
    A_str = A.to_csv(sep="\t", header=False, index=False)
    B_str = B.to_csv(sep="\t", header=False, index=False)

    A_bed = pybedtools.BedTool(A_str, from_string=True)
    B_bed = pybedtools.BedTool(B_str, from_string=True)

    # Intersect: returns TP CRIPSRi peaks that overlap with CLR peaks
    overlap =  A_bed.intersect(B_bed, wa=True)  

    # Convert result to DataFrame
    overlap_df = overlap.to_dataframe(names=["chr", "start", "end"])
    overlap_df[["start", "end"]] = overlap_df[["start", "end"]].astype(int)
    
    return overlap_df
    

# Find the peaks (Chr, Start, End, CLR_value) associated with gene #
def EG_pair_by_name(gene, CLR_Matrix):
    
    if gene not in CLR_Matrix.columns:
        print(f"Warning: {gene} not found in CLR matrix.")
    else:
        EG_pair = CLR_Matrix[gene].copy()
        EG_pair = EG_pair[EG_pair > 0]
        
        
        GATA1_EG_pair = EG_pair.index.str.extract(r'^(chr\w+):(\d+)-(\d+)$')
        GATA1_EG_pair.columns = ['chr', 'start', 'end']
        GATA1_EG_pair['CLR_value'] = EG_pair.values


        return GATA1_EG_pair


# ---------------------------------------------CRRSPRi Functions------------------------------------------------



'''


# convert final CLR cCRE matrixes into list 
def convert_columns_to_list_if_not_zero_sum(df):
    # Dictionary to store columns with non-zero sum as lists
    column_lists = []
    
    # Iterate over each column
    for column in df.columns:
        col_sum = df[column].sum()
        
        # Check if the sum is not zero
        if col_sum != 0:
            column_lists.append(column)
    
    return column_lists


# convert final CLR cCRE matrix into list 
def convert_columns_to_list_chr(CLR_final):
    column_lists = {}
    
    for key, value in CLR_final.items():
        list = convert_columns_to_list_if_not_zero_sum(value)
        column_lists[key] = list

    return column_lists




# find overlaps between two df file ( but in df format)
def count_overlapping_peaks(test_df, reference_df):
    
    # Convert test_df and reference_df into pybedtools BedTool objects
    test_bed = pybedtools.BedTool.from_dataframe(test_df)
    reference_bed = pybedtools.BedTool.from_dataframe(reference_df)
    
    # Perform intersection (overlap) between the test and reference intervals
    overlap = reference_bed.intersect(test_bed, u=True)
    
    # Count the number of overlaps
    overlap_count = len(overlap)
    
    return overlap, overlap_count



# find overlaps between a given bed file and dict containing chromosome CREs
# overlap: overlapped sections
# counts: each chr, the number of overlapped 
# chr_counts: orginal number of chr section
# total_counts: adding all counts
# total_chr_count: adding all chr_count
def count_overlapping_chr(test_dict, bed_file_path):
    bed = pd.read_csv(bed_file_path, sep = '\t', comment = '#')
    bed = bed.iloc[:, [0, 1, 2, 9]]
    bed.columns = ['Chr', 'Start', 'End','Dnase']
    new_bed = bed[bed['Dnase']=='DNase-only']
    
    overlaps = {}
    counts = {}
    chr_count = {}
    total_chr_count = 0
    total_counts = 0
    
    for key,  value in test_dict.items():
        bed_chr = new_bed[new_bed['Chr'] == key]
        
        df = pd.DataFrame(value, columns=['Coordinates'])

        # Split the 'Coordinates' column into 'Chromosome', 'Start', 'End'
        df[['Chromosome', 'Position']] = df['Coordinates'].str.split(':', expand=True)
        df[['Start', 'End']] = df['Position'].str.split('-', expand=True)

        # Drop the original 'Coordinates' and 'Position' columns (optional)
        df = df.drop(columns=['Coordinates', 'Position'])
        
        overlap, count = count_overlapping_peaks(bed_chr,df)
        
        overlaps[key] = overlap
        counts[key] = count
        chr_count[key] = len(df.index)
        
        total_counts = total_counts + count
        total_chr_count = total_chr_count + len(df.index)
         
        
    
    return overlaps,counts, chr_count, total_counts, total_chr_count


'''