version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_pgen2bed.wdl" as pgen_conversion
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks

workflow basic_admixture {
  input {
    Array[File] vcf
    File? pop
    Int n_ancestral_populations
    Boolean cross_validation = false
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
				vcf = file
		}

		if (prune_variants) {
			call variant_tasks.pruneVars {
				input:
					pgen = subsetVariants.subset_pgen,
					pvar = subsetVariants.subset_pvar,
					psam = subsetVariants.subset_psam,
					window_size = window_size,
					shift_size = shift_size,
					r2_threshold = r2_threshold
			}
		}

		File subset_pgen = select_first([pruneVars.out_pgen, subsetVariants.subset_pgen])
		File subset_pvar = select_first([pruneVars.out_pvar, subsetVariants.subset_pvar])
		File subset_psam = select_first([pruneVars.out_psam, subsetVariants.subset_psam])
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				pgen = subset_pgen,
				pvar = subset_pvar,
				psam = subset_psam
		}
	}

	File merged_pgen = select_first([mergeFiles.out_pgen, pruneVars.out_pgen[0], subsetVariants.subset_pgen[0]])
	File merged_pvar = select_first([mergeFiles.out_pvar, pruneVars.out_pvar[0], subsetVariants.subset_pvar[0]])
	File merged_psam = select_first([mergeFiles.out_psam, pruneVars.out_psam[0], subsetVariants.subset_psam[0]])

  	if (remove_relateds) {
		call sample_tasks.removeRelateds {
			input:
				pgen = merged_pgen,
				pvar = merged_pvar,
				psam = merged_psam,
				max_kinship_coefficient = max_kinship_coefficient
		}
	}

	File final_pgen = select_first([removeRelateds.out_pgen, merged_pgen])
	File final_pvar = select_first([removeRelateds.out_pvar, merged_pvar])
	File final_psam = select_first([removeRelateds.out_psam, merged_psam])

  call pgen_conversion.pgen2bed {
    input:
      pgen = final_pgen,
      pvar = final_pvar,
      psam = final_psam
  }

  if (defined(pop)) {
    call subset_pop {
      input:
        fam = pgen2bed.out_fam,
        pop = select_first([pop])
    }
  }
  File? this_pop = if (defined(pop)) then subset_pop.out_pop else pop

  call Admixture_t {
    input:
      bed = pgen2bed.out_bed,
      bim = pgen2bed.out_bim,
      fam = pgen2bed.out_fam,
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
