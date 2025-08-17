# =================================================================================================
#     Dependencies
# =================================================================================================


configfile: "config.yaml"


## this rule runs a python script that will generate the bed file containing the contributing loci and their effect sizes based on values of dn alpha
## the bed file will be then used to annotate a vcf file that will be used by SliM to run the simulations


rule all:
    input:
        expand(
            "results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_ecotype_counts.csv",
            pi=config["pi"],
            outcrossing_rate=config["outcrossing_rate"],
            selection=config["selection"],
            heritability=config["heritability"],
            replicates_arq=config["replicates_arq"],
            optima=config["optima"],
            replicates_sim=config["replicates_sim"],    
        ),


rule build_population_for_sim:
    input:
        og_tree_offset=config["og_tree_offset"],
        og_vcf_offset=config["og_vcf_offset"],
    output:
        tree_seq_causalloci="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/tree_seq_causalloci.trees",
        loci_effectsize="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/loci_effectsize.csv",
        phenotypes="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/phenotypes.csv",
    params:
        pi=lambda wildcards: str(wildcards.pi),
        replicates_arq=lambda wildcards: str(wildcards.replicates_arq),
        beta=config["beta"],
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/outxing_{outcrossing_rate}_arq_pi{pi}_{replicates_arq}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/build_population_for_sim.py"

rule run_slim_simulation:
    input:
        tree_seq_causalloci="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/tree_seq_causalloci.trees",
    output: 
        output_tree_gen6=temp("results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen6.trees"),
        #output_tree_gen10="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen10.trees",
        output_pop_size_early="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_pop_size_early.txt",
        output_va="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_va.txt",
        output_mfitness="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_mfitness.txt",
        output_vfitness="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vfitness.txt",
        output_mpheno="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_mpheno.txt",
        output_vpheno="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vpheno.txt",
        #output_new_optimum="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_new_optimum.txt",
        #output_adj_variance="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_adj_variance.txt",
        output_maxphenotype="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_maxphenotype.txt",
        output_minphenotype="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_minphenotype.txt",

    resources:
        mem_mb=40960,
    conda:
        "envs/base_env.yaml"
    benchmark:
        "benchmarks/outxing_{outcrossing_rate}_arq_pi{pi}_{replicates_arq}_{heritability}_{selection}_optima{optima}_subp{replicates_sim}.txt"
    script:
        "scripts/slim.sh"

rule tree_postprocessing:
    input:
        og_tree_offset=config["og_tree_offset"],
        mapper_ids=config['mapper_realid_metadataid'],
        output_sim_tree="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_tree_output_gen6.trees",
    output:
        output_vcf=temp("results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcf_gen6.vcf"),
    resources:
        mem_mb=30720,
        limit_space=1,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/tree_postprocessing.py"

    rule calc_ecotype_counts:
    input:
        nonhet_pos=config['nonhet_pos'],
        og_vcf_offset=config["og_vcf_offset"],
        ecotypes_grenenet=config['ecotypes_grenenet'],
        output_vcf_offset="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_vcf_gen6.vcf",
    output:
        ecotype_counts="results/outxing_{outcrossing_rate}/arq_pi{pi}_{replicates_arq}/{heritability}/{selection}/optima{optima}/subp{replicates_sim}_ecotype_counts.csv",
    resources:
        mem_mb=30720,
    threads: 20,
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/calc_ecotype_counts.py"
