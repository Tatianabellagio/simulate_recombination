import os

pheno = snakemake.input['pheno_file'] 

output_gen6_trees = pheno.replace('_minphenotype.txt', '_tree_output_gen6.trees')


def create_empty_file_if_not_exists(filename):
    if not os.path.exists(filename):
        # Create an empty file
        with open(filename, 'w') as file:
            pass

create_empty_file_if_not_exists(output_gen6_trees)