"""
Tree Files Filling Script
=========================

This script creates empty tree files when populations die early in simulations.
It ensures the pipeline doesn't break due to missing tree sequence files.

Purpose:
- Handle cases where SLiM simulations fail or populations go extinct
- Create placeholder files so downstream analysis can continue
- Prevent pipeline failures due to missing tree sequence outputs
"""

import os

# Get input and output files from Snakemake
pheno_file = snakemake.input['pheno_file']
sim_treegen3 = snakemake.output['sim_treegen3']
sim_treegen10 = snakemake.output['sim_treegen10']

def create_empty_file_if_not_exists(filename):
    """
    Create an empty file if it doesn't exist.
    
    Args:
        filename: Path to the file to create
    """
    try:
        if not os.path.exists(filename):
            # Create an empty file
            with open(filename, 'w') as file:
                pass
            print(f"Created empty file: {filename}")
        else:
            print(f"File already exists: {filename}")
    except Exception as e:
        print(f"Error creating file {filename}: {e}")

# Create empty tree files for both generation 3 and 10
# This handles cases where populations die early and don't produce these files
print(f"Processing phenotype file: {pheno_file}")
print("Creating empty tree files if they don't exist...")

create_empty_file_if_not_exists(sim_treegen3)
create_empty_file_if_not_exists(sim_treegen10)

print("Tree files filling completed successfully!")