version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/variant_filtering.wdl" as variant_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/sample_filtering.wdl" as sample_tasks
import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/main/file_tasks.wdl" as file_tasks

workflow prep_target {
	input {
		File ref_bim
		Array[File] vcf
		Boolean remove_relateds = true
		Int kinship_degree_filter = 3
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
		call sample_tasks.king_ibdseg {
			input: 
				bed = merged_bed,
				bim = merged_bim,
				fam = merged_fam,
				degree = kinship_degree_filter
		}

		call sample_tasks.findRelated {
			input: 
				king_file = king_ibdseg.kin0, 
				estimator = "PropIBD", 
				degree = kinship_degree_filter
		}

		if (findRelated.has_relatives) {
			call removeSamples {
				input: 
				bed = merged_bed,
				bim = merged_bim,
				fam = merged_fam,
				samples_to_remove = findRelated.related_samples,
				suffix = "unrel"
			}
		}
	}

	File final_bed = select_first([removeSamples.out_bed, merged_bed])
	File final_bim = select_first([removeSamples.out_bim, merged_bim])
	File final_fam = select_first([removeSamples.out_fam, merged_fam])

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

task removeSamples {
    input {
        File bed
        File bim
        File fam
        File samples_to_remove
        String suffix = "subset"
        Int mem_gb = 16
    }

    Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 10
    String basename = basename(bed, ".bed")

    command <<<
        command="plink2 --bed ~{bed} --bim ~{bim} --fam ~{fam} \
        --remove ~{samples_to_remove} \
        --output-chr chrM \
        --make-bed \
        --out ~{basename}_~{suffix}"
        printf "${command}\n"
        ${command}
    >>>

    output {
        File out_bed="~{basename}_~{suffix}.bed"
        File out_bim="~{basename}_~{suffix}.bim"
        File out_fam="~{basename}_~{suffix}.fam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}
