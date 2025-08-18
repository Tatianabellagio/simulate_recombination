import os

output_gen6_trees = snakemake.output['sim_tree_gen6'] 
print(output_gen6_trees)
print(os.path.exists(output_gen6_trees))

def create_empty_file_if_not_exists(output_gen6_trees):
    if not os.path.exists(output_gen6_trees):
        # Create an empty file
        print(f"Creating empty file: {output_gen6_trees}")
        with open(output_gen6_trees, 'w') as file:
            pass

create_empty_file_if_not_exists(output_gen6_trees)