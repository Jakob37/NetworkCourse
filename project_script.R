#!/usr/bin/env Rscript

message("Loading packages")
library(igraph)
library(tidyverse)
library(gridExtra)

outdir <- "figs/"
dir.create(outdir)
message("Creating output dir: ", outdir)

corr_fp <- "p_long_p0001.tsv.gz"
message("Loading transcript correlations: ", corr_fp)
trans_corr <- read_tsv(corr_fp)
# trans_corr <- read_tsv("p_long_p0000001.tsv.gz")

prot_fp <- "ppi_human_stringdb_ensp2ensg2.tab.gz"
message("Loading protein interactions: ", prot_fp)
prot_int <- read_tsv(prot_fp)

message("Transcript entries: ", nrow(trans_corr))
message("Prot entries: ", nrow(prot_int))
message("Transcripts per prot: ", nrow(trans_corr) / nrow(prot_int))

message("Generating graphs")

#thresholds <- c(1e-3, 1e-4, 1e-5)
thresholds <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12)

prot_graph <- graph_from_data_frame(prot_int[, c("gene1", "gene2")])
prot_deg <- igraph::degree(prot_graph, mode="out")

plts <- list()
corrs <- list()

message("Processing correlations for thresholds: ", paste(thresholds, collapse=", "))
for (i in seq_len(length(thresholds))) {

    thres <- thresholds[i]
    message("Processing threshold: ", thres)
    
    subset_trans <- trans_corr %>% filter(value < thres)
    trans_graph <- graph_from_data_frame(subset_trans)
    trans_deg <- igraph::degree(trans_graph, mode="out")
    common_names <- intersect(names(trans_deg), names(prot_deg))

    message("Common entries found: ", length(common_names))

    common_trans_out <- trans_deg[common_names]
    common_prot_out <- prot_deg[common_names]
    out_deg_df <- data.frame(
        trans=common_trans_out,
        prot=common_prot_out
    )
    
    plts[[i]] <- grid.arrange(
        ggplot(out_deg_df, aes(x=trans)) + geom_histogram(bins=100) + theme_classic() + xlim(NA, 2000) + ggtitle(paste0("Common transcripts, outdegrees")),
        ggplot(out_deg_df, aes(x=prot)) + geom_histogram(bins=100) + theme_classic() + xlim(NA, 2000) + ggtitle(paste0("Common proteins, outdegrees"))
    )

    message("Pearson corr for thres ", thres, ": ", cor(out_deg_df$trans, out_deg_df$prot, method="pearson"))
    message("Spearman corr for thres ", thres, ": ", cor(out_deg_df$trans, out_deg_df$prot, method="spearman"))
    
    print(cor.test(out_deg_df$trans, out_deg_df$prot, method="pearson"))
    print(cor.test(out_deg_df$trans, out_deg_df$prot, method="spearman"))
    
}

message("Saving figures to: ", outdir)
for (i in seq_len(length(thresholds))) {
    thres <- thresholds[i]
    ggsave(plts[[i]], file=paste0(outdir, "/hists_", thres, ".png"))
}
























