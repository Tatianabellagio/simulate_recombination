"""
Ecotype Counts Analysis Script
==============================

This script analyzes simulated VCF files to count how many individuals belong to each ecotype
in the simulated populations. It compares the simulated genotypes against the original 
GreneNet ecotype classifications to track population structure changes.

Purpose:
- Track how different mating systems (outcrossing rates) affect population structure
- Monitor genetic diversity changes over generations  
- Compare simulated populations to the original GreneNet diversity
- Generate data for downstream analysis of selection effects on ecotype distribution

Input files:
- output_vcf_offset: VCF from SLiM simulations (after 6 generations)
- nonhet_pos: Positions that are non-heterozygous in the original data
- og_vcf_offset: Original GreneNet VCF file (reference)
- ecotypes_grenenet: Mapping of individuals to ecotypes

Output:
- ecotype_counts: CSV file with counts of each ecotype in the simulated population
"""

import numpy as np
import pandas as pd
import allel
import os
from collections import defaultdict
import multiprocessing

## Get input/output files from Snakemake
## This script processes VCFs from a particular combination of architecture, heritability, and selection
## and expands across all replicates and optima
output_vcf_offset = snakemake.input['output_vcf_offset']  # Simulated VCF from SLiM
nonhet_pos = snakemake.input['nonhet_pos']                # Non-heterozygous positions
og_vcf_offset = snakemake.input['og_vcf_offset']          # Original GreneNet VCF
ecotypes_grenenet = snakemake.input['ecotypes_grenenet']  # Ecotype mapping
ecotype_counts = snakemake.output['ecotype_counts']        # Output file

def filtering_pos(nonhet_pos, pos_new, geno_og, geno_new):
    """
    Filter genotypes to only include positions that are non-heterozygous in the original data.
    
    This ensures we're comparing genotypes at the same genetic positions between
    the original GreneNet data and the simulated populations.
    
    Args:
        nonhet_pos: Array of non-heterozygous positions from original data
        pos_new: Positions in the new (simulated) VCF
        geno_og: Genotype array from original VCF
        geno_new: Genotype array from new VCF
    
    Returns:
        geno_og_rpos: Filtered original genotypes at common positions
        geno_new_rpos: Filtered new genotypes at common positions
    """
    # Find positions that exist in both datasets
    pos_to_keep = np.intersect1d(nonhet_pos, pos_new)
    
    # Create masks to filter genotypes at these common positions
    mask_pos_ogvcf = pd.Series(pos_og).isin(pos_to_keep)
    geno_og_rpos = geno_og[mask_pos_ogvcf]
    
    mask_pos_newvcf = pd.Series(pos_new).isin(pos_to_keep)
    geno_new_rpos = geno_new[mask_pos_newvcf]
    
    return geno_og_rpos, geno_new_rpos

def get_ecotype_geno_mapper(geno_og_rpos):
    """
    Create a mapping from genotype patterns to ecotype names.
    
    This function creates a dictionary where each unique genotype pattern
    (as bytes) maps to its corresponding ecotype from the original GreneNet data.
    
    Args:
        geno_og_rpos: Filtered genotype array from original VCF
    
    Returns:
        ecotype_geno_mapper: Dictionary mapping genotype bytes to ecotype names
    """
    # Transpose genotypes so each row represents one individual
    geno_og_rpos = np.swapaxes(geno_og_rpos, 0, 1)
    
    ecotype_geno_mapper = {}
    # For each individual, create a genotype pattern and map it to their ecotype
    for i, j in zip(geno_og_rpos, samples):
        geno = i.tobytes()  # Convert genotype array to bytes for hashing
        ecotype_geno_mapper[geno] = j
    
    return ecotype_geno_mapper

def get_ecotype_counts(geno_new_rpos, pop_name):
    """
    Count how many individuals in the simulated population belong to each ecotype.
    
    This function analyzes the simulated genotypes and maps them back to
    the original ecotype classifications to track population structure changes.
    
    Args:
        geno_new_rpos: Filtered genotype array from simulated VCF
        pop_name: Name identifier for the population (for output column)
    
    Returns:
        ecotype_countsdf: DataFrame with ecotype counts for this population
    """
    # Transpose genotypes so each row represents one individual
    geno_new_rpos = np.swapaxes(geno_new_rpos, 0, 1)
    
    # Initialize a defaultdict to store the genotype counts
    ecotype_counts = defaultdict(int)
    
    # Count genotypes in the simulated population
    for i in geno_new_rpos:
        sample = i.tobytes()  # Convert genotype to bytes for lookup
        # Look up the ecotype for this genotype, default to 'other' if not found
        ecotype = ecotype_geno_mapper.get(sample, 'other')
        ecotype_counts[ecotype] += 1
    
    # Create DataFrame with ecotype counts
    name = pop_name
    ecotype_countsdf = pd.DataFrame(list(ecotype_counts.items()), columns=['ecotype', name])
    
    # Clean up ecotype names (remove suffixes after underscore)
    ecotype_countsdf['ecotype'] = ecotype_countsdf['ecotype'].str.split('_').str[0]
    
    return ecotype_countsdf

## Load the non-heterozygous positions data
nonhet_pos = np.array(pd.read_csv(nonhet_pos))

## Load the original GreneNet VCF file
vcf_og = allel.read_vcf(og_vcf_offset, fields=['calldata/GT', 'variants/POS', 'samples'])
geno_og = vcf_og['calldata/GT']      # Genotype data
pos_og = vcf_og['variants/POS']      # Position data
samples = vcf_og['samples']          # Sample names (ecotypes)

## Note: The following code was commented out but shows the original approach
## of loading ecotype mappings directly from a CSV file
#ecotypes_grenenet = pd.read_csv(ecotypes_grenenet, dtype=object)
#ecotypes_grenenet.columns= ['ecotype']
#ecotypes_grenenet = pd.concat([ecotypes_grenenet, pd.DataFrame(data = {'ecotype': ['other']}, index=[231])],axis=0)

print(f"Processing VCF: {output_vcf_offset}")

## Extract population name from the file path for identification
name = output_vcf_offset.split('/')[-2] + '_' + output_vcf_offset.split('/')[-1][0:5]

## Check if the VCF file exists and has content
if os.path.exists(output_vcf_offset) and os.path.getsize(output_vcf_offset) <= 1:
    # File exists but is empty (likely simulation failed)
    pass
elif os.path.exists(output_vcf_offset) and os.path.getsize(output_vcf_offset) > 1:
    print('Processing non-empty VCF file')
    
    ## Import the simulated VCF file
    vcf_new = allel.read_vcf(output_vcf_offset, fields=['calldata/GT', 'variants/POS'])
    
    ## Extract positions and genotype arrays from simulated VCF
    pos_new = vcf_new['variants/POS']      # Variant positions
    geno_new = vcf_new['calldata/GT']      # Genotype data
    
    ## Filter genotypes to common positions between original and simulated data
    geno_og_rpos, geno_new_rpos = filtering_pos(nonhet_pos, pos_new, geno_og, geno_new)
    
    ## Create the genotype-to-ecotype mapping based on filtered positions
    ecotype_geno_mapper = get_ecotype_geno_mapper(geno_og_rpos)
    
    ## Count ecotypes in the simulated population
    ecotype_countsdf = get_ecotype_counts(geno_new_rpos, name)
    print(f"Ecotype counts for {name}:")
    print(ecotype_countsdf)
    
    ## Note: The following code was commented out but shows how to merge
    ## results with previous analyses
    #ecotypes_grenenet = ecotypes_grenenet.merge(ecotype_countsdf, how='left', on ='ecotype')
    #print(ecotypes_grenenet)

## Save the results or create empty file if no data
try:
    print(f"Saving ecotype counts to: {ecotype_counts}")
    ecotype_countsdf.to_csv(ecotype_counts)
except NameError:
    print('No ecotype data to save - creating empty file')
    # Create an empty file with the same name if ecotype_countsdf is not defined
    # This happens when simulations fail or produce empty VCFs
    open(ecotype_counts, 'w').close()