version 1.0 

workflow custom_plotting {
    input {
        File ancestry_frac
        File cluster_groups
    }

    call plot {
        input: 
            ancestry_frac = ancestry_frac, 
            cluster_groups = cluster_groups
    }

    output {
        File ancestry_plot = plot.ancestry_plot
    }
}

task plot {
    input {
        File ancestry_frac
        File cluster_groups
    }

    command <<<
        R << RSCRIPT
        library(tidyverse)
        library(RColorBrewer)

        ancestry_frac <- read_table("~{ancestry_frac}", col_names = FALSE)
        K <- ncol(ancestry_frac) - 1
        names(ancestry_frac) <- c('sample_id', paste0('K', 1:K))

        cluster_map <- read_tsv("~{cluster_groups}", col_names = c('new', 'old'))

        print(head(ancestry_frac))
        print(head(cluster_map))

        ancestry_frac <- ancestry_frac %>% rename(!!!setNames(cluster_map[['new']], cluster_map[['old']]))

        cluster_order <- paste0("K", 1:K)
        ancestry_frac <- ancestry_frac %>% arrange(across(all_of(cluster_order)))
        ancestry_frac <- ancestry_frac %>% mutate(n = row_number())

        ancestry_frac <- ancestry_frac %>% 
            pivot_longer(
                cols = starts_with("K"),
                names_to = "Cluster",
                values_to = "K"
            )

        ancestry_frac$Cluster <- factor(ancestry_frac$Cluster, levels = cluster_order)

        d2 <- brewer.pal(8, 'Set1')
        s2 <- brewer.pal(8, 'Pastel1')
        colormap <- setNames(c(d2, s2)[1:K], paste0("K", 1:K))

        p <- ggplot(ancestry_frac, aes(x=n, y=K, fill=Cluster, color=Cluster)) + 
            geom_bar(stat='identity') + 
                scale_fill_manual(values=colormap, breaks=rev(names(colormap))) + 
                scale_color_manual(values=colormap, breaks=rev(names(colormap))) + 
                theme_classic() + 
                theme(
                    axis.line=element_blank(), 
                    axis.ticks.x=element_blank(), 
                    axis.text.x=element_blank(), 
                    axis.title.x=element_blank(), 
                    axis.ticks.y=element_blank(), 
                    axis.text.y=element_blank(), 
                    axis.title.y=element_blank(), 
                    panel.spacing=unit(0, 'in')
                )
        ggsave("ancestry_plot.png", p, width=10, height=4)
        RSCRIPT
    >>>

    output {
        File ancestry_plot = "ancestry_plot.png"
    }

    runtime {
        docker: "rocker/tidyverse:4.3.1"  # R + tidyverse
        memory: "8G"
        cpu: 1
    }
}