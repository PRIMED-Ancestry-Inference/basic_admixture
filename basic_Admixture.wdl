version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks

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

	call plot_admixture {
		input:
			ancestry_frac = Admixture_t.ancestry_fractions
	}

	output {
		File ancestry_fractions = Admixture_t.ancestry_fractions
		File allele_frequencies = Admixture_t.allele_frequencies
		File ancestry_plot = plot_admixture.plot
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
		writeLines(dat[['pop']], '~{outfile}'); \
		"
	>>>

	output {
		File out_pop = outfile
	}

	runtime {
		docker: "rocker/tidyverse:4"
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


task plot_admixture {
	input {
		File ancestry_frac
	}

	command <<<
		Rscript -e "\
		library(tidyverse); \
		library(RColorBrewer); \
		dat <- read_delim('~{ancestry_frac}', col_names=FALSE); \
		K <- ncol(dat) - 1; \
		names(dat) <- c('sample_id', paste0('K', 1:K)); \
		dat <- arrange(dat, across(starts_with('K'))); \
		dat <- mutate(dat, n=row_number()); \
		dat <- pivot_longer(dat, starts_with('K'), names_to='Cluster', values_to='K'); \
		d2 <- brewer.pal(8, 'Dark2'); s2 <- brewer.pal(8, 'Set2'); \
		colormap <- setNames(c(d2, s2)[1:K], paste0('K', 1:K)); \
		ggbar <- ggplot(dat, aes(x=n, y=K, fill=Cluster, color=Cluster)) + \
		geom_bar(stat='identity') + \
		scale_fill_manual(values=colormap, breaks=rev(names(colormap))) + \
		scale_color_manual(values=colormap, breaks=rev(names(colormap))) + \
		theme_classic() + \
		theme(axis.line=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), panel.spacing=unit(0, 'in')); \
		ggsave('admixture_plot.png', width=11, height=4); \
		"
	>>>

	output {
		File plot = "admixture_plot.png"
	}

	runtime {
		docker: "rocker/tidyverse:4"
	}
}

