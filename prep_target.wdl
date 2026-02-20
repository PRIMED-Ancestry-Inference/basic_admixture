version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks

workflow prep_target {
	input {
		File ref_bim
		Array[File] vcf
		Boolean remove_relateds = true
		Float? max_kinship_coefficient
		Int mem_gb = 16
	}

	call selectColumn {
		input:
			ref_variants = ref_bim,
			variant_id_col = 2
	}

	scatter (file in vcf) {
		call variant_tasks.subsetVariants {
			input:
				vcf = file,
				variant_files = [selectColumn.id_file],
				output_chr = "26"
		}
	}

	if (length(vcf) > 1) {
		call file_tasks.mergeFiles {
			input:
				bed = subsetVariants.subset_bed,
				bim = subsetVariants.subset_bim,
				fam = subsetVariants.subset_fam,
				output_chr = "26"
		}
	}

	File merged_bed = select_first([mergeFiles.out_bed, subsetVariants.subset_bed[0]])
	File merged_bim = select_first([mergeFiles.out_bim, subsetVariants.subset_bim[0]])
	File merged_fam = select_first([mergeFiles.out_fam, subsetVariants.subset_fam[0]])

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

	output {
		File bed = final_bed
		File bim = final_bim
		File fam = final_fam
	}
}


task selectColumn {
	input {
		File ref_variants
		Int variant_id_col = 2
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
