import os

output_gen6_trees = snakemake.output['sim_tree_gen6'] 

def create_empty_file_if_not_exists(filename):
    if not os.path.exists(filename):
        # Create an empty file
        print(f"Creating empty file: {filename}")
        with open(filename, 'w') as file:
            pass

create_empty_file_if_not_exists(output_gen6_trees)