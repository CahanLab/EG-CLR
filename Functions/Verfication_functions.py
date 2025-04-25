import pandas as pd
from muon import atac as ac
from muon import prot as pt
import anndata as ad
import random
random.seed(10)
import pybedtools


### verfication


# ---------------------------------------------CRRSPRi Functions------------------------------------------------


# compare CLR prediction with CRISPRi data #
def CRISPRi_comparison(CRISPRi,CLR):
    
    # Find TP for CRISPRi
    CRISPRi_TP = CRISPRi[CRISPRi['Significant'] == True].copy()
    CRISPRi_TP = CRISPRi_TP[['chr', 'start', 'end']].copy()
    
    CLR_str = CLR.to_csv(sep="\t", header=False, index=False)
    CRISPRi_TP_str = CRISPRi_TP.to_csv(sep="\t", header=False, index=False)

    CLR_bed = pybedtools.BedTool(CLR_str, from_string=True)
    CRISPRi_TP_bed = pybedtools.BedTool(CRISPRi_TP_str, from_string=True)

    # Intersect: returns TP CRIPSRi peaks that overlap with CLR peaks
    overlap =  CRISPRi_TP_bed.intersect(CLR_bed, wa=True)  

    # Convert result to DataFrame
    overlap_df = overlap.to_dataframe(names=["chr", "start", "end"])
    
    return overlap_df
    

# Find detected peaks in the CRISPRi data #
def find_detected_peaks(CLR, CRISPRi): 
    
    # limit to chr, start and end
    CLR = CLR[["chr", "start", "end"]].copy()
    CRISPRi_new = CRISPRi[["chr", "start", "end"]].copy()

    # Load as BedTool
    CLR_str = CLR.to_csv(sep="\t", header=False, index=False)
    CRISPRi_str = CRISPRi_new.to_csv(sep="\t", header=False, index=False)

    CLR_bed = pybedtools.BedTool(CLR_str, from_string=True)
    CRISPRi_bed = pybedtools.BedTool(CRISPRi_str, from_string=True)
    
    # Drop non-detected peaks in CRISPRi
    CRISPRi_filtered_bed = CRISPRi_bed.intersect(CLR_bed, wa=True)
    detected_peaks = CRISPRi_filtered_bed.to_dataframe(names=["chr", "start", "end"])
    
    # Add the original columns
    df_merged = pd.merge(detected_peaks,CRISPRi, on=['chr', 'start', 'end'], how='left')

    return df_merged


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


# find overlaps between two bed file ( but in df format)
def count_overlapping_peaks(test_df, reference_df):
    
    # Convert test_df and reference_df into pybedtools BedTool objects
    test_bed = pybedtools.BedTool.from_dataframe(test_df)
    reference_bed = pybedtools.BedTool.from_dataframe(reference_df)
    
    # Perform intersection (overlap) between the test and reference intervals
    overlap = reference_bed.intersect(test_bed, u=True)
    
    # Count the number of overlaps
    overlap_count = len(overlap)
    
    return overlap, overlap_count


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


'''
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