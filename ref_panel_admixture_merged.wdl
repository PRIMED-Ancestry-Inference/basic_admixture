version 1.0 

import "https://raw.githubusercontent.com/UW-GAC/primed-bcftools/5243ee37ea5e360788a9ba43fe24cb83a27292bf/extract_vcf_ids.wdl" as extract_vcf_ids
import "basic_admixture_merged.wdl" as basic_admixture
import "projected_admixture_merged.wdl" as projected_admixture

workflow ref_panel_admixture_merged {
    input {
        Array[File] study_vcf_file
        Array[File] ref_vcf_file
        Int n_ancestral_populations
        File? ref_sample
        File ref_pop
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
            sample_file = ref_sample,
            pop = ref_pop
        # output = ancestry_fractions, allele_frequencies, ancestry_plot, bed, bim, fam 
    }

    call projected_admixture.projected_admixture {
        input:
            ref_allele_freq = basic_admixture.allele_frequencies,
		    vcf = study_vcf_file
        # output: bed, bim, fam 
    }

    call merge {
        input: 
            ref_bed = basic_admixture.bed,
            ref_bim = basic_admixture.bim,
            ref_fam = basic_admixture.fam,

            proj_bed = projected_admixture.bed,
            proj_bim = projected_admixture.bim,
            proj_fam = projected_admixture.fam
        # output: merged_fam, merged_vcf 
    }

    call make_merged_pop_file {
        input: 
            merged_fam = merge.merged_fam,
            ref_ancestry_frac = basic_admixture.ancestry_fractions, 
            ref_pop = ref_pop
        # output sample_file, pop_file
    }

    call basic_admixture.basic_admixture as merged_admixture {
        input:
            vcf = [merge.merged_vcf],
            n_ancestral_populations = n_ancestral_populations,
            pop = make_merged_pop_file.pop_file,
            sample_file = make_merged_pop_file.sample_file
    }

    call plot_admixture {
        input: 
            ancestry_frac = merged_admixture.ancestry_fractions,
            proj_fam = projected_admixture.fam,
            ref_pop = ref_pop
    }
    output {
        File ancestry_fractions = merged_admixture.ancestry_fractions
		File allele_frequencies = merged_admixture.allele_frequencies
		File plot = plot_admixture.plot
    }
}

task merge {
    input {
        File ref_bed
        File ref_bim
        File ref_fam

        File proj_bed
        File proj_bim
        File proj_fam

        Int mem_gb = 16
    }

    Int disk_size = ceil(
        2.5 * (
            size(ref_bed, "GB") +
            size(ref_bim, "GB") +
            size(ref_fam, "GB") +
            size(proj_bed, "GB") +
            size(proj_bim, "GB") +
            size(proj_fam, "GB")
        )
    ) + 20

    command <<<
        set -e -o pipefail

        plink \
        --bed ~{ref_bed} --bim ~{ref_bim} --fam ~{ref_fam} \
        --bmerge ~{proj_bed} ~{proj_bim} ~{proj_fam} \
        --make-bed \
        --out tmp

        plink \
        --bfile tmp \
        --recode vcf-iid bgz \
        --out merged_combined
    >>>

    output {
        File merged_vcf = "merged_combined.vcf.gz"
        File merged_fam = "tmp.fam"
    }

    runtime {
        docker: "quay.io/biocontainers/plink:1.90b6.21--h516909a_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}

task make_merged_pop_file {
  input {
    File merged_fam
    File ref_ancestry_frac
    File ref_pop
  }

  command <<<
        R << RSCRIPT
        library(tidyverse)

        fam <- read.table('~{merged_fam}', stringsAsFactors=F)
        sample <- data.frame(FID = fam[[1]], IID = fam[[2]])
        sample_tmp <- data.frame(IID = fam[,2])

        pop_tmp <- read.table('~{ref_pop}', stringsAsFactors=F, col.names = c('IID', 'POP'))

        #left join sample and pop, fill in blanks with "-"
        pop <- left_join(sample_tmp, pop_tmp, by = 'IID') %>% mutate(POP = if_else(is.na(POP), '-', POP))

        write.table(sample, 'sample_file.txt', quote=F, row.names=F, col.names=F)
        write.table(pop, 'pop_file.txt', quote=F, row.names=F, col.names=F)
        RSCRIPT
    >>>

  output {
    File sample_file = "sample_file.txt"
    File pop_file = "pop_file.txt"
  }

  runtime {
    docker: "rocker/tidyverse:4"
  }
}

task plot_admixture {
    input {
        File ancestry_frac
        File proj_fam
        File ref_pop
    }

	command <<<
        Rscript -e "\
        library(tidyverse); \
        library(RColorBrewer); \

        fam <- read_table('~{proj_fam}', col_names=FALSE); \
        target_ids <- fam[[2]]; \

        pop <- read_table('~{ref_pop}', col_names=c('sample_id','POP')); \

        dat <- read_delim('~{ancestry_frac}', col_names=FALSE); \
        K <- ncol(dat) - 1; \
        names(dat) <- c('sample_id', paste0('K', 1:K)); \

        dat_full <- left_join(dat, pop, by='sample_id'); \
        cluster_means <- dat_full %>% \
            filter(POP != '-') %>% \
            group_by(POP) %>% \
            summarise(across(starts_with('K'), mean)); \
        write.table(cluster_means, 'cluster_means.txt', quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t'); \

        dat <- dat %>% filter(sample_id %in% target_ids); \
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
        File cluster_means = "cluster_means.txt"
	}

	runtime {
		docker: "rocker/tidyverse:4"
	}
}