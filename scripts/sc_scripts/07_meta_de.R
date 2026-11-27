# Compares this dataset's pseudobulk DESeq2 results (05_tdp43_de_pseudobulk.R)
# with the bulk TDP-43- vs TDP-43+ nuclei results (06_liu_de.R).

suppressMessages({
  library(tidyverse)
  library(ggplot2)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

# NEEDS INPUT: matches the cell type tested in the 05_tdp43_de_* scripts.
cell_type <- "L2/3 EN"
cell_type_file <- str_replace_all(cell_type, "[ /]", "_")

results_dir <- "results/sc_results/07_meta_de/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

snrna <- read.csv(paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/deseq2.csv"))
bulk <- read.csv("results/sc_results/06_liu_de/deseq.csv")

bulk <- bulk[, c(1, 3, 7)]
snrna <- snrna[, c(1, 3, 7)]

# for MAST results only:
# snrna <- snrna %>% 
#   mutate(avg_log2FC = avg_log2FC * -1)

colnames(bulk) <- c("gene", "bulk_fc", "bulk_padj")
colnames(snrna) <- c("gene", "snrna_fc", "snrna_padj")

df <- full_join(bulk, snrna,
                by = "gene") %>% 
  na.omit()

head(df)

df <- df %>%
  mutate(bulk_deg = if_else(bulk_padj < 0.05 & abs(bulk_fc) > 0.25, T, F),
         snrna_deg = if_else(snrna_padj < 0.05 & abs(snrna_fc) > 0.25, T, F)) %>%
  mutate(deg = case_when(bulk_deg == T & snrna_deg == T ~ "Both",
                         bulk_deg == T & snrna_deg == F ~ "Bulk",
                         bulk_deg == F & snrna_deg == T ~ "snRNA",
                         .default = "Not DEG"))

p <- df %>%
  filter(deg != "Not DEG") %>%
  ggplot(aes(x = bulk_fc,
             y = snrna_fc)) + 
  geom_point(aes(color = deg)) + 
  theme_linedraw()
ggsave(p,
       filename = paste0(results_dir, cell_type_file, "_vs_bulk_lfc.png"),
       units = "in", dpi = 600,
       height = 6, width = 7)
