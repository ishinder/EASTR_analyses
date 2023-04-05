# plot reference-matching vs novel junctions using a junctions table

library(tidyverse)
library(ggplot2)
library(httpgd)
library(extrafont)
font_import()
font_import(pattern="arial.ttf", prompt=FALSE)

group_col <- "Aligner"

# source data
junction_file <- "/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra/eastr_juncs_alns_summary.tsv"
#junction_file <- "/ccb/salz2/shinder/projects/EASTR_tests2/arabidopsis_methylation/eastr_juncs_alns_summary.tsv"

# parse
junction_stats <- read_tsv(junction_file) %>%
    mutate(`EASTR Reference Junctions` = `Reference Junctions` - `Removed Reference Junctions`) %>%
    mutate(`EASTR Novel Junctions` = `Novel Junctions` - `Removed Novel Junctions`)
groups <- junction_stats[[group_col]] %>% unique()


plot_junctions <- function(groups) {
    print(paste0(group_col, ": ", groups))

    group_stats <- junction_stats[junction_stats[[group_col]] %in% groups,]
    subset_samples <- group_stats[["SRR ID"]][1:3]
    subset_stats <- group_stats[group_stats[["SRR ID"]] %in% subset_samples,]

    print(
        ggplot(data=subset_stats) +
            geom_point(aes(x="Novel Junctions", y="Reference Junctions", color=Aligner)) +
            geom_point(aes(x="EASTR Novel Junctions", y="EASTR Reference Junctions", color=Aligner)) + 
            facet_wrap(vars(`SRR ID`)) + theme_bw() +
            theme(text = element_text(family = "Roboto Mono"))
        )
}


hgd()
hgd_browse()

plot_junctions(groups)

dev.off()
