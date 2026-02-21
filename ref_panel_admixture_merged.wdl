version 1.0 

import "https://raw.githubusercontent.com/UW-GAC/primed-bcftools/5243ee37ea5e360788a9ba43fe24cb83a27292bf/extract_vcf_ids.wdl" as extract_vcf_ids
import "prep_ref.wdl" as prep_ref
import "prep_target.wdl" as prep_target

workflow ref_panel_admixture_merged {
    input {
        Array[File] study_vcf_file
        Array[File] ref_vcf_file
        Int n_ancestral_populations
        File? ref_sample
        File ref_pop
        Boolean cross_validation = false
    }

    call extract_vcf_ids.extract_vcf_ids {
        input:
            vcf_file = study_vcf_file
        # output = variant_file
    }

    call prep_ref.prep_ref {
        input:
            vcf = ref_vcf_file,
            ref_variants = extract_vcf_ids.variant_file,
            n_ancestral_populations = n_ancestral_populations,
            sample_file = ref_sample,
            pop = ref_pop
        # output = bed, bim, fam, ref_pop
    }

    call prep_target.prep_target {
        input:
            ref_bim = prep_ref.bim,
		    vcf = study_vcf_file
        # output: bed, bim, fam 
    }

    call merge {
        input: 
            ref_bed = prep_ref.bed,
            ref_bim = prep_ref.bim,
            ref_fam = prep_ref.fam,

            proj_bed = prep_target.bed,
            proj_bim = prep_target.bim,
            proj_fam = prep_target.fam
        # output: merged_bed, merged_bim, merged_fam 
    }

    call make_merged_pop_file {
        input: 
            merged_fam = merge.merged_fam,
            ref_pop = prep_ref.ref_pop
        # output sample_file, pop_file
    }

    call Admixture_t {
		input:
			bed = merge.merged_bed,
			bim = merge.merged_bim,
			fam = merge.merged_fam,
			pop = make_merged_pop_file.pop_file,
			n_ancestral_populations = n_ancestral_populations,
			cross_validation = cross_validation
	}

    call plot_admixture {
        input: 
            ancestry_frac = Admixture_t.ancestry_fractions,
            proj_fam = prep_target.fam,
            ref_pop = prep_ref.ref_pop
    }
    output {
        File ancestry_fractions = Admixture_t.ancestry_fractions
		File allele_frequencies = Admixture_t.allele_frequencies
		File plot = plot_admixture.plot
        File cluster_means = plot_admixture.cluster_means
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

        #plink \
        #--bfile tmp \
        #--recode vcf-iid bgz \
        #--out merged_combined
    >>>

    output {
        #File merged_vcf = "merged_combined.vcf.gz"
        File merged_bed = "tmp.bed"
        File merged_bim = "tmp.bim"
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

task Admixture_t {
    input {
        File bed
        File bim
        File fam
        File? pop # two column, ID and pop
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
        if [ -f ~{pop} ]; then ln -s ~{pop} ~{basename}.pop.full; fi
        cut -f2 ~{basename}.pop.full > ~{basename}.pop
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