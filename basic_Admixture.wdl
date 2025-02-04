version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/filter_chr/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/filter_chr/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/filter_chr/file_tasks.wdl" as file_tasks

workflow basic_admixture {
	input {
		Array[File] vcf
		File? ref_variants
		File? sample_file
		File? pop
		Int n_ancestral_populations
		Boolean cross_validation = false
		Int? genome_build
		Boolean prune_variants = true
		Boolean remove_relateds = true
		Float? min_maf
		Float? max_kinship_coefficient
		Int? window_size
		Int? shift_size
		Int? r2_threshold
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_files = select_all([ref_variants]),
				sample_file = sample_file,
				genome_build = genome_build,
				min_maf = min_maf
		}

		if (prune_variants) {
			call variant_tasks.pruneVars {
				input:
					bed = subsetVariants.subset_bed,
					bim = subsetVariants.subset_bim,
					fam = subsetVariants.subset_fam,
					window_size = window_size,
					shift_size = shift_size,
					r2_threshold = r2_threshold,
					output_chr = "26"
			}
		}

		File subset_bed = select_first([pruneVars.out_bed, subsetVariants.subset_bed])
		File subset_bim = select_first([pruneVars.out_bim, subsetVariants.subset_bim])
		File subset_fam = select_first([pruneVars.out_fam, subsetVariants.subset_fam])
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				bed = subset_bed,
				bim = subset_bim,
				fam = subset_fam,
				output_chr = "26"
		}
	}

	File merged_bed = select_first([mergeFiles.out_bed, pruneVars.out_bed[0], subsetVariants.subset_bed[0]])
	File merged_bim = select_first([mergeFiles.out_bim, pruneVars.out_bim[0], subsetVariants.subset_bim[0]])
	File merged_fam = select_first([mergeFiles.out_fam, pruneVars.out_fam[0], subsetVariants.subset_fam[0]])

	if (remove_relateds) {
		call sample_tasks.removeRelateds {
			input:
				bed = merged_bed,
				bim = merged_bim,
				fam = merged_fam,
				max_kinship_coefficient = max_kinship_coefficient,
				output_chr = "26"
		}
	}

	File final_bed = select_first([removeRelateds.out_bed, merged_bed])
	File final_bim = select_first([removeRelateds.out_bim, merged_bim])
	File final_fam = select_first([removeRelateds.out_fam, merged_fam])

	if (defined(pop)) {
		call subset_pop {
			input:
				fam = final_fam,
				pop = select_first([pop])
		}
	}
	File? this_pop = if (defined(pop)) then subset_pop.out_pop else pop

	call Admixture_t {
		input:
			bed = final_bed,
			bim = final_bim,
			fam = final_fam,
			pop = this_pop,
			n_ancestral_populations = n_ancestral_populations,
			cross_validation = cross_validation
	}

	output {
		File ancestry_fractions = Admixture_t.ancestry_fractions
		File allele_frequencies = Admixture_t.allele_frequencies
	}
}


task subset_pop {
	input{
		File fam
		File pop # two columns: ID pop
	}

	String outfile = basename(fam) + ".pop"

	command <<<
		Rscript -e "\
		library(readr); \
		library(dplyr); \
		fam <- read_delim('~{fam}', col_types='-c----', col_names='id'); \
		dat <- read_delim('~{pop}', col_names=c('id', 'pop')); \
		dat <- left_join(fam, dat); \
		dat <- mutate(dat, pop=ifelse(is.na(pop), '-', pop)); \
		writeLines(dat[['pop']], '~{outfile}'); \
		"
	>>>

	output {
		File out_pop = outfile
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.16.0"
	}
}


task Admixture_t {
	input {
		File bed
		File bim
		File fam
		File? pop
		File? P # include this for use with projected_admixture
		Int n_ancestral_populations
		Boolean cross_validation = false
		Int mem_gb = 16
		Int n_cpus = 4
	}

	Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
	String basename = basename(bed, ".bed")

	command <<<
		set -e -o pipefail
		ln -s ~{bed} ~{basename}.bed
		ln -s ~{bim} ~{basename}.bim
		ln -s ~{fam} ~{basename}.fam
		if [ -f ~{pop} ]; then ln -s ~{pop} ~{basename}.pop; fi
		if [ -f ~{P} ]; then ln -s ~{P} ~{basename}.~{n_ancestral_populations}.P.in; fi
		/admixture_linux-1.3.0/admixture ~{if defined(P) then "-P" else ""} ~{if cross_validation then "--cv" else ""} \
			~{basename}.bed ~{n_ancestral_populations} ~{if defined(pop) then "--supervised" else ""} \
			-j~{n_cpus}
		paste -d' ' <(cut -f2 ~{basename}.fam) ~{basename}.~{n_ancestral_populations}.Q > ~{basename}.~{n_ancestral_populations}.ancestry_frac
		paste -d' ' <(cut -f2 ~{basename}.bim) ~{basename}.~{n_ancestral_populations}.P > ~{basename}.~{n_ancestral_populations}.allele_freq
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/admixture_docker:v1.0.0"
		disks: "local-disk " + disk_size + " SSD"
		memory: mem_gb + " GB"
		cpu: n_cpus
	}

	output {
		File ancestry_fractions = "~{basename}.~{n_ancestral_populations}.ancestry_frac"
		File allele_frequencies = "~{basename}.~{n_ancestral_populations}.allele_freq"
	}
}
