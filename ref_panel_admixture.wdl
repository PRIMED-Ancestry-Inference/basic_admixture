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
            n_ancestral_populations = n_ancestral_populations

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