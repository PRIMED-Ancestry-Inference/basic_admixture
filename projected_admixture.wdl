version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_pgen2bed.wdl" as pgen_conversion
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks
import "basic_Admixture.wdl" as admixture

workflow projected_admixture {
	input{
		File ref_allele_freq
    	Array[File] vcf
    	Boolean cross_validation = false
	}

	call selectColumn {
		input:
			ref_variants = ref_allele_freq,
			variant_id_col = 1
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = selectColumn.id_file
		}
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				pgen = subsetVariants.subset_pgen,
				pvar = subsetVariants.subset_pvar,
				psam = subsetVariants.subset_psam
		}
	}

	File final_pgen = select_first([mergeFiles.out_pgen, subsetVariants.subset_pgen[0]])
	File final_pvar = select_first([mergeFiles.out_pvar, subsetVariants.subset_pvar[0]])
	File final_psam = select_first([mergeFiles.out_psam, subsetVariants.subset_psam[0]])

	call pgen_conversion.pgen2bed {
		input:
			pgen = final_pgen,
			pvar = final_pvar,
			psam = final_psam
	}

	call admixReady {
		input:
			ref_allele_freq = ref_allele_freq,
			bim = pgen2bed.out_bim
	}
	
	call admixture.Admixture_t {
		input:
			bed = pgen2bed.out_bed,
			bim = pgen2bed.out_bim,
			fam = pgen2bed.out_fam,
			P = admixReady.subset_P,
			n_ancestral_populations = admixReady.k,
			cross_validation = cross_validation
	}

	output {
		File ancestry_fractions = Admixture_t.ancestry_fractions
		File allele_frequencies = Admixture_t.allele_frequencies
	}

	meta {
    author: "Jonathan Shortt"
    email: "jonathan.shortt@cuanschutz.edu"
    description: "This workflow is used to project a genetic test dataset (in VCF format) into clusters (\"ancestral populations\") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected into the clusters determined by the .P."
	}
}


task selectColumn {
	input {
		File ref_variants
		Int variant_id_col = 1
	}

	command <<<
		cut -f ~{variant_id_col} ~{ref_variants} > variant_ids.txt
	>>>

	output {
		File id_file = "variant_ids.txt"
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.17.0"
	}
}


task admixReady {
	input{
		File ref_allele_freq
		File bim
	}
	
	String basename = basename(ref_allele_freq)

	command <<<
		#subset ref_allele_freq to only variants also in bim, print all but first column to get .P file
		awk 'FNR==NR{a[$2]; next}{if($1 in a){print $0}}' ~{bim} ~{ref_allele_freq} | cut -d' ' -f2- > ~{basename}_admixReady.P

		my_k=$(head -n 1 ~{basename}_admixReady.P | awk '{print NF}')
		printf "\nProjection will be run with k=${my_k} clusters.\n"
		printf "${my_k}" > tmp.txt

		snp_count=$(wc -l < ~{basename}_admixReady.P)
		printf "\nProjection uses ${snp_count} snps.\n"
	>>>
	
	output {
		File subset_P = "~{basename}_admixReady.P"
		String k = read_string("tmp.txt")
	}
	
	runtime {
    	docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	}
}
