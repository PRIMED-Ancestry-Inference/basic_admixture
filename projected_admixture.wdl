version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks
import "basic_Admixture.wdl" as admixture

workflow projected_admixture {
	input {
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
				variant_files = [selectColumn.id_file]
		}
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				bed = subsetVariants.subset_bed,
				bim = subsetVariants.subset_bim,
				fam = subsetVariants.subset_fam
		}
	}

	File final_bed = select_first([mergeFiles.out_bed, subsetVariants.subset_bed[0]])
	File final_bim = select_first([mergeFiles.out_bim, subsetVariants.subset_bim[0]])
	File final_fam = select_first([mergeFiles.out_fam, subsetVariants.subset_fam[0]])

	call admixReady {
		input:
			ref_allele_freq = ref_allele_freq,
			bim = final_bim
	}

	call admixture.Admixture_t {
		input:
			bed = final_bed,
			bim = final_bim,
			fam = final_fam,
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
		set -e -o pipefail
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
