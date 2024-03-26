# basic_admixture

Workflows to run [ADMIXTURE](https://dalexander.github.io/admixture/) Borrows some code from [here](github.com/broadinstitute/palantir-workflows/Admixture).

## basic_admixture

This workflow is used to run ADMIXTURE on a genetic dataset (in VCF format). Optionally,
the data is pruned for linkage equilibrium, and related individuals are removed.

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
pop | Optional file with known population labels. This should be a two-column file in the format "id pop"
n_ancestral_populations | number of clusters to infer
cross_validation | Boolean for whether run run cross-validation (default false)
prune_variants | Boolean for whether to do LD pruning on the variants (default true)
remove_relateds | Boolean for whether to remove samples with relatedness above max_kinship_coefficient (default true)
max_kinship_coefficient | if remove_relateds is true, remove one of each pair of samples with kinship > this value (default 0.0442 = 3rd degree relatives)
window_size | window size for LD pruning (default 10,000)
shift_size | shift size for LD pruning (default 1000)
r2_threshold | r2 threshold for LD pruning (default 0.1)


Outputs:

output | description
--- | ---
ancestry_fractions | ancestry fractions for each sample (.Q output from ADMIXTURE with sample IDs as first column)
allele_frequencies | allele frequencies for each variants (.P output from ADMIXTURE with variant IDs as first column)


## projected_admixture

This workflow is used to project a genetic test dataset (in VCF format) into clusters ("ancestral populations") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected into the clusters determined by the .P.
	

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
ref_allele_freq | allele_frequency output from basic_admixture workflow (first column variant id, remaining colums are frequencies in each cluster)
cross_validation | Boolean for whether run run cross-validation (default false)


Outputs:

output | description
--- | ---
ancestry_fractions | ancestry fractions for each sample (.Q output from ADMIXTURE with sample IDs as first column)
allele_frequencies | allele frequencies for each variants (.P output from ADMIXTURE with variant IDs as first column)
