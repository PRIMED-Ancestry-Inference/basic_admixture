version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/pgen2bed/plink2_pgen2bed.wdl" as pgen_conversion
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/vcf_input/create_pca_projection.wdl" as tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/vcf_input/projected_pca.wdl" as file_tasks

workflow projected_admixture {
	input{
		File P
		File ref_variants
    	Array[File] vcf
		Int? seed # https://wdl-docs.readthedocs.io/en/latest/WDL/different_parameters/
		Boolean? cv
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
	
	call run_admixture_projected {
		input:
			bed = pgen2bed.out_bed,
			bim = pgen2bed.out_bim,
			fam = pgen2bed.out_fam,
			P = admixReady.subset_P,
			k = summary.k,
			seed = seed,
			cv = cv
	}
	
	meta {
    author: "Jonathan Shortt"
    email: "jonathan.shortt@cuanschutz.edu"
    description: "This workflow is used to project a genetic test dataset (in plink format, i.e., .bed/.bim/.fam) into clusters (\"ancestral populations\") using ADMIXTURE. First, the cluster file (.P produced by ADMIXTURE) and the test dataset are both subset to contain the same set of variants (Note: this workflow assumes that variants from both the .P and test dataset have been previously harmonized such that variants follow the same naming convention, alleles at each site are ordered identically, and variants are sorted). Then the test dataset is projected into the clusters determined by the .P."
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

task run_admixture_projected {
	input {
    	File bed
    	File bim
    	File fam
    	File P
    	Int k
    	Int? seed
    	Boolean cv=false
    	Int mem_gb = 16 #https://github.com/openwdl/wdl/pull/464
    	Int n_cpus = 4
  	}

    Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")
	#ln --symbolic ${P} ${basename}.${k}.P.in

	command <<<
		ln --symbolic ~{P} ~{basename}.~{k}.P.in
		command='/admixture_linux-1.3.0/admixture ~{if (cv) then "--cv" else ""} ~{if defined(seed) then "-s ~{seed}" else "-s time"} -j~{n_cpus} -P ~{bed} ~{k}'
		printf "${command}\n"
		${command} | tee ~{basename}_projection.~{k}.log
		#/admixture_linux-1.3.0/admixture ~{if (cv) then "--cv" else ""} ~{if defined(seed) then "-s ~{seed}" else "-s time"} -j~{n_cpus} -P ~{bed} ~{k} | tee ~{basename}_projection.~{k}.log
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/admixture_docker:v1.0.0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
  	}

	output {
		File ancestry_fractions = "~{basename}.~{k}.Q"
		File allele_frequencies = "~{basename}.~{k}.P"
		File admixture_log = "~{basename}_projection.~{k}.log"
	}
}
