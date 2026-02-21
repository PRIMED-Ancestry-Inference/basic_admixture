version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks

workflow prep_ref {
	input {
		Array[File] vcf
		File? ref_variants
		File? sample_file
		File pop
		Int n_ancestral_populations
		Int? genome_build
		Boolean prune_variants = true
		Boolean remove_relateds = true
		Float? min_maf
		Float? max_kinship_coefficient
		Int? window_size
		Int? shift_size
		Float? r2_threshold
	}

	if (defined(ref_variants)) {
		call remove_chr_prefix {
			input: 
				variant_file = select_first([ref_variants, ""])
		}
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_files = select_all([remove_chr_prefix.output_file]),
				sample_file = sample_file,
				genome_build = genome_build,
				min_maf = min_maf,
				output_chr = "26"
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

	call subset_pop {
		input:
			fam = final_fam,
			pop = select_first([pop])
	}

	output {
		File bed = final_bed
		File bim = final_bim
		File fam = final_fam
		File ref_pop = subset_pop.out_pop
	}
}


task remove_chr_prefix {
	input {
		File variant_file
	}

	command <<<
		sed 's/^chr//' ~{variant_file} > chr_int.txt
	>>>

	output {
		File output_file = "chr_int.txt"
	}

	runtime {
		docker: "rocker/tidyverse:4"
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
		write.table(dat, 'out_pop.txt', quote=FALSE, sep = '\t', col.names=F, row.names=FALSE) \
		"
	>>>

	output {
		File out_pop = "out_pop.txt"
	}

	runtime {
		docker: "rocker/tidyverse:4"
	}
}


