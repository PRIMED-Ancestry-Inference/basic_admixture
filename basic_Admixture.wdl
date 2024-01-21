version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/plink2_vcf2pgen.wdl" as vcf_conversion
import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/pgen2bed/plink2_pgen2bed.wdl" as pgen_conversion
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/vcf_input/create_pca_projection.wdl" as subset_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/vcf_input/projected_pca.wdl" as file_tasks

workflow basic_admixture {
  input {
    Array[File] vcf
    File? pop
    Int n_ancestral_populations
    Boolean cross_validation = false
    Boolean prune_variants = true
    Boolean remove_relateds = true
    Float? max_kinship_coefficient
		Int? window_size
		Int? shift_size
		Int? r2_threshold
  }

	scatter (file in vcf) {
		call vcf_conversion.vcf2pgen {
			input:
				vcf_file = file
		}

		if (prune_variants) {
			call subset_tasks.pruneVars {
				input:
					pgen = vcf2pgen.out_pgen,
					pvar = vcf2pgen.out_pvar,
					psam = vcf2pgen.out_psam,
					window_size = window_size,
					shift_size = shift_size,
					r2_threshold = r2_threshold
			}
		}

		File subset_pgen = select_first([pruneVars.out_pgen, vcf2pgen.out_pgen])
		File subset_pvar = select_first([pruneVars.out_pvar, vcf2pgen.out_pvar])
		File subset_psam = select_first([pruneVars.out_psam, vcf2pgen.out_psam])
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				pgen = subset_pgen,
				pvar = subset_pvar,
				psam = subset_psam
		}
	}

	File merged_pgen = select_first([mergeFiles.out_pgen, pruneVars.out_pgen[0], vcf2pgen.out_pgen[0]])
	File merged_pvar = select_first([mergeFiles.out_pvar, pruneVars.out_pvar[0], vcf2pgen.out_pvar[0]])
	File merged_psam = select_first([mergeFiles.out_psam, pruneVars.out_psam[0], vcf2pgen.out_psam[0]])

  	if (remove_relateds) {
		call subset_tasks.removeRelateds {
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

  call Admixture_t {
    input:
      bed = pgen2bed.out_bed,
      bim = pgen2bed.out_bim,
      fam = pgen2bed.out_fam,
      pop = pop,
      n_ancestral_populations = n_ancestral_populations,
      cross_validation = cross_validation
  }

  output {
    File ancestry_fractions = Admixture_t.ancestry_fractions
    File allele_frequencies = Admixture_t.allele_frequencies
    File reference_variants = pgen2bed.out_bim
  }
}


task Admixture_t {
  input {
    File bed
    File bim
    File fam
    File? pop
    Int n_ancestral_populations
    Boolean cross_validation = false
    Int mem = 16
    Int n_cpus = 4
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
  String basename = basename(bed, ".bed")

  command <<<
    ln -s ~{bed} ~{basename}.bed
    ln -s ~{bim} ~{basename}.bim
    ln -s ~{fam} ~{basename}.fam
    if [ -f ~{pop} ]; then ln -s ~{pop} ~{basename}.pop; fi
    /admixture_linux-1.3.0/admixture ~{if cross_validation then "--cv" else ""} \
      ~{basename}.bed ~{n_ancestral_populations} ~{if defined(pop) then "--supervised" else ""} \
      -j~{n_cpus}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/admixture_docker:v1.0.0"
		disks: "local-disk " + disk_size + " SSD"
    memory: mem + " GB"
    cpu: n_cpus
  }

  output {
    File ancestry_fractions = "~{basename}.~{n_ancestral_populations}.Q"
    File allele_frequencies = "~{basename}.~{n_ancestral_populations}.P"
  }
}
