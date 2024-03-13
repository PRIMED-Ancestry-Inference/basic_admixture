version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_pgen2bed.wdl" as pgen_conversion
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/create_pca_projection.wdl" as tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/projected_pca.wdl" as file_tasks
import "basic_Admixture.wdl" as admixture

workflow projected_admixture {
	input{
		File P
		File ref_variants
    	Array[File] vcf
    	Boolean cross_validation = false
	}
	
	call tasks.identifyColumns {
		input:
			ref_variants = ref_variants
	}

	scatter (file in vcf) {
		call file_tasks.subsetVariants {
			input:
				vcf = file,
				variant_file = ref_variants,
				variant_id_col = identifyColumns.id_col
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
			P = P, 
			ref_variants = ref_variants,
			variant_id_col = identifyColumns.id_col,
			pvar = final_pvar
	}
	
	call summary{
		input:
			P = P,
			snps = admixReady.snps
	}
	
	call admixture.Admixture_t {
		input:
			bed = pgen2bed.out_bed,
			bim = pgen2bed.out_bim,
			fam = pgen2bed.out_fam,
			P = admixReady.subset_P,
			n_ancestral_populations = summary.k,
			cross_validation = cross_validation
	}
	
	meta {
    author: "Jonathan Shortt"
    email: "jonathan.shortt@cuanschutz.edu"
    description: "This workflow is used to project a genetic test dataset (in VCF format) into clusters (\"ancestral populations\") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected into the clusters determined by the .P."
	}

  output {
    File ancestry_fractions = Admixture_t.ancestry_fractions
    File allele_frequencies = Admixture_t.allele_frequencies
  }
}

task admixReady {
	input{
		File P
		File ref_variants
		Int variant_id_col
		File pvar
	}
	
	String basename = basename(P, ".P")

	command <<<
		#get a list of variant names, save to extract.txt
		cut -f3 ~{pvar} > extract.txt

		#extract in-common variants from ref.P
		paste -d'\t' <(cut -f~{variant_id_col} ~{ref_variants}) ~{P} > tmp
		head tmp
		head ~{ref_variants}
		head ~{P}
		awk 'FNR==NR{a[$1]; next}{if($1 in a){print $0}}' extract.txt tmp | cut -f~{variant_id_col}- > ~{basename}_admixReady.P
	>>>
	
	output {
		File subset_P="~{basename}_admixReady.P"
		File snps="extract.txt"
	}
	
	runtime {
    	docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	}
}

task summary {
	input{
		File P
		File snps
	}
	
	#String my_k = ceil(awk NR==1{print NF} P) #maybe do a check on the file to make sure all lines are the same to make more robust?
	
	command <<<
		my_k=$(head -n 1 ~{P} | awk '{print NF}')
		printf "\nProjection will be run with k=${my_k} clusters.\n"
		printf "${my_k}" > tmp.txt
		snp_count=$(wc -l < ~{snps})
		printf "\nProjection uses ${snp_count} snps.\n"
	>>>
	
	output {
		String k = read_string("tmp.txt")
	}
	
	runtime{
		docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
	}
}
