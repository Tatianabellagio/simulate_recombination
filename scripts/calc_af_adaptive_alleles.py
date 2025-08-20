import numpy as np
import pandas as pd
import allel

POS_CSV = snakemake.input['loci_effectsize']  # Simulated VCF from SLiM
VCF_PATH = snakemake.input['vcf']                # Non-heterozygous positions
af_recomb = snakemake.output['af_recomb']         # Original GreneNet VCF

HET_THRESH = 30  # threshold to label recombinant vs non-recombinant

# 1) read targets (adaptive positions) and VCF (only what we need)
pos = pd.read_csv(POS_CSV)["pos"].to_numpy()
callset = allel.read_vcf(VCF_PATH, fields=["variants/POS", "calldata/GT", "samples"])

variants = callset["variants/POS"]
gt = allel.GenotypeArray(callset["calldata/GT"])
samples = np.array(callset.get("samples", [f"S{i}" for i in range(gt.shape[1])]))

# 2) basic genotype-derived arrays
n_alt = gt.to_n_alt()                  # (n_variants, n_samples) with 0/1/2
is_het = gt.is_het()
het_per_sample = is_het.sum(axis=0)    # genome-wide heterozygous count per individual

# 3) select adaptive variants
mask = np.isin(variants, pos)
idx_adapt = np.flatnonzero(mask)

# 4) split samples by heterozygosity threshold
idx_R  = np.flatnonzero(het_per_sample >  HET_THRESH)   # recombinant
idx_NR = np.flatnonzero(het_per_sample <= HET_THRESH)   # non-recombinant

# 5) ALT allele copies contributed by each group (vector per variant)
alt_R  = n_alt[np.ix_(idx_adapt, idx_R)].sum(axis=1)
alt_NR = n_alt[np.ix_(idx_adapt, idx_NR)].sum(axis=1)
alt_tot = alt_R + alt_NR

# 6) frequencies and shares
ploidy = 2
N_R, N_NR = len(idx_R), len(idx_NR)
N_tot = N_R + N_NR

AF_tot = alt_tot / (ploidy * N_tot) if N_tot else np.zeros_like(alt_tot, dtype=float)
AF_R   = alt_R   / (ploidy * N_R)   if N_R   else np.zeros_like(alt_R,   dtype=float)
AF_NR  = alt_NR  / (ploidy * N_NR)  if N_NR  else np.zeros_like(alt_NR,  dtype=float)

share_R  = np.divide(alt_R,  alt_tot, out=np.zeros_like(alt_R,  dtype=float), where=alt_tot>0)
share_NR = np.divide(alt_NR, alt_tot, out=np.zeros_like(alt_NR, dtype=float), where=alt_tot>0)

# 7) tidy output
out = pd.DataFrame({
    "pos": variants[idx_adapt],
    "AF_total": AF_tot,
    "AF_recombinant": AF_R,
    "AF_nonrecombinant": AF_NR,
    "share_alt_from_recombinant": share_R,
    "share_alt_from_nonrecombinant": share_NR,
    "alt_copies_total": alt_tot,
    "alt_copies_recombinant": alt_R,
    "alt_copies_nonrecombinant": alt_NR,
    "N_recombinant": N_R,
    "N_nonrecombinant": N_NR,
})

out.to_csv(af_recomb)