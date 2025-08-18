# Recombination effects on adaptive evolution in Arabidopsis

We are investigating the effects of recombination through outcrossing on adaptive evolution in Arabidopsis
### Two Experimental Conditions:
1. **With Recombination**: Outcrossing rates of 0%, 5%, 10%, and 15% (based on literature)
2. **Without Recombination**: 0% outcrossing (purely clonal) as control

### Simulation Parameters:
- **Selection Strength**: .1
- **Heritability**: 0.7
- **Environmental Change**: 1,2,3 std
- **Population Size**: cap at 900 at the end of each generation 
- **Generations**: 5 generations to observe evolutionary change

## Repository Structure

```
.
├── Snakefile              # Main Snakemake workflow for simulations
├── config.yaml           # Configuration: outcrossing rates, selection, heritability
├── scripts/              # Core simulation and analysis scripts
│   ├── build_population_for_sim.py    # Creates initial populations from GreneNet data
│   ├── slim.sh                        # Runs SLiM simulations with parameters
│   ├── arabidopsis_evolve_treeseq.slim # SLiM script with outcrossing logic
│   ├── tree_postprocessing.py         # Converts tree sequences to VCF
│   └── calc_ecotype_counts.py        # Calcualtes the number of different ecotypes and outcrossers based on initial vcf 
├── treeseq/             # Tree sequence conversion scripts
├── analysis/            # Analysis notebooks and scripts
└── results/             # Simulation outputs organized by parameters
```