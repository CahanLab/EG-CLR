import pybedtools


### verfication

# find overlaps between two bed file ( but in df format)
def count_overlapping_peaks(test_df, reference_df):
    """
    Counts how many peaks from reference_df overlap with test_df using pybedtools.
    
    Parameters:
    test_df (pd.DataFrame): Test DataFrame with columns ['chr', 'start', 'end']
    reference_df (pd.DataFrame): Reference DataFrame with columns ['chr', 'start', 'end']
    
    Returns:
    int: The number of overlapping peaks from reference_df with test_df
    """
    
    # Convert test_df and reference_df into pybedtools BedTool objects
    test_bed = pybedtools.BedTool.from_dataframe(test_df)
    reference_bed = pybedtools.BedTool.from_dataframe(reference_df)
    
    # Perform intersection (overlap) between the test and reference intervals
    overlap = reference_bed.intersect(test_bed, u=True)
    
    # Count the number of overlaps
    overlap_count = len(overlap)
    
    return overlap, overlap_count

