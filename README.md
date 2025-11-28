# basic_admixture

Workflows to run [ADMIXTURE](https://dalexander.github.io/admixture/) Borrows some code from [here](github.com/broadinstitute/palantir-workflows/Admixture).

## basic_admixture

This workflow is used to run ADMIXTURE on a genetic dataset (in VCF format). Either unsupervised or supervised mode is supported.
Optionally, the data is pruned for linkage equilibrium, and related individuals are removed.

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
ref_variants | Optional file with variants to use.
sample_file | Optional file with sample IDs to keep. This should be in the form accepted by PLINK ("FID\tIID"), so the sample ID column should be duplicated.
pop | Optional file with known population labels. This should be a two-column file in the format "id pop". If this file is provided, ADMIXTURE is run in supervised mode; otherwise, the clustering is unsupervised.
n_ancestral_populations | number of clusters to infer
cross_validation | Boolean for whether to run cross-validation (default false)
prune_variants | Boolean for whether to do LD pruning on the variants (default true)
min_maf | minimum MAF for variants to include (optional)
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

This workflow is used to project a genetic test dataset (in VCF format) into clusters ("ancestral populations") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted.) Then the test dataset is projected into the clusters determined by the .P.
	

Inputs:

input | description
--- | ---
vcf | Array of VCF files (possibly split by chromosome)
ref_allele_freq | allele_frequency output from basic_admixture workflow (first column variant id, remaining colums are frequencies in each cluster)
cross_validation | Boolean for whether to run cross-validation (default false)


Outputs:

output | description
--- | ---
ancestry_fractions | ancestry fractions for each sample (.Q output from ADMIXTURE with sample IDs as first column)
allele_frequencies | allele frequencies for each variants (.P output from ADMIXTURE with variant IDs as first column)



## ref_panel_admixture

This workflow runs extract_vcf_ids, basic_admixture and projected_admixture in succesion. This workflow is used to extract variant identifiers from a set of study VCFs (extract_vcf_ids), runs ADMIXTURE using the reference panel VCFs and the variants identified in the first step (basic_admixture), and projects the study VCFs (projected_admixture)
	

Inputs:

input | description
--- | ---
study_vcf_file | Array of study VCF files (possibly split by chromosome)
ref_vcf_file | Array of reference panel VCF files (possibly split by chromosome)
n_ancestral_populations | number of clusters to infer


Outputs:

output | description
--- | ---
ref_ancestry_fractions | ancestry fractions for each sample in reference panel (.Q output from ADMIXTURE with sample IDs as first column)
ref_allele_frequencies | allele frequencies for each variants in reference panel (.P output from ADMIXTURE with variant IDs as first column)
ref_ancestry_plot | ancestry plot for samples in reference panel 
ancestry_fractions | ancestry fractions for each sample in study data (.Q output from ADMIXTURE with sample IDs as first column)
allele_frequencies | allele frequencies for each variants in study data (.P output from ADMIXTURE with variant IDs as first column)
ancestry_plot | ancestry plot for samples in study data



## custom_plotting

This workflow makes ancestry plots with cluster labels from input file. 
	

Inputs:

input | description
--- | ---
ancestry_frac | File with ancestry fractions 
cluster_groups | File with two columns. First column is new cluster labels, the second column is which original cluster label maps to the new cluster label


Outputs:

output | description
ancestry_plot | ancestry plot for samples with new labels 
