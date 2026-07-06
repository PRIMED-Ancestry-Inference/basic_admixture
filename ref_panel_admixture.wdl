version 1.0 

import "https://raw.githubusercontent.com/UW-GAC/primed-bcftools/5243ee37ea5e360788a9ba43fe24cb83a27292bf/extract_vcf_ids.wdl" as extract_vcf_ids
import "basic_Admixture.wdl" as basic_admixture
import "projected_admixture.wdl" as projected_admixture

workflow ref_panel_admixture{
    input {
        Array[File] study_vcf_file
        Array[File] ref_vcf_file
        Int n_ancestral_populations
        Int mem_gb = 16
        File? pop
        Boolean cross_validation = false
        Int? genome_build
        Boolean prune_variants = true
        Boolean remove_relateds = true
        Float? min_maf
        Float? max_kinship_coefficient
        Int? window_size
        Int? shift_size
        Float? r2_threshold
    }

    call extract_vcf_ids.extract_vcf_ids {
        input:
            vcf_file = study_vcf_file

        # output = variant_file
    }

    call basic_admixture.basic_admixture {
        input:
            vcf = ref_vcf_file,
            ref_variants = extract_vcf_ids.variant_file,
            n_ancestral_populations = n_ancestral_populations,
            pop = pop,
            cross_validation = cross_validation,
            genome_build = genome_build,
            prune_variants = prune_variants,
            remove_relateds = remove_relateds,
            min_maf = min_maf,
            max_kinship_coefficient = max_kinship_coefficient,
            window_size = window_size,
            shift_size = shift_size,
            r2_threshold = r2_threshold
        # output = ancestry_fractions, allele_frequencies, ancestry_plot
    }

    call projected_admixture.projected_admixture {
        input:
            ref_allele_freq = basic_admixture.allele_frequencies,
		    vcf = study_vcf_file, 
            mem_gb = mem_gb

        # output = ancestry_fractions, allele_frequencies, ancestry_plot
    }

    output { 
        File ref_ancestry_fractions = basic_admixture.ancestry_fractions
        File ref_allele_frequencies = basic_admixture.allele_frequencies
        File ref_ancestry_plot = basic_admixture.ancestry_plot
        File ancestry_fractions = projected_admixture.ancestry_fractions
        File allele_frequencies = projected_admixture.allele_frequencies
        File ancestry_plot = projected_admixture.ancestry_plot
    }
}
