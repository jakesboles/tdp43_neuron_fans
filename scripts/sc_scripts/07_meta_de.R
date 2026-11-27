library(tidyverse)
library(ggplot2)

setwd("/projects/b1169/boles/JSB149")

snrna <- read.csv("tab_data/05_tdp43_de/deseq2.csv")
bulk <- read.csv("tab_data/06_liu_de/deseq.csv")

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

df %>%
  filter(deg != "Not DEG") %>%
  ggplot(aes(x = bulk_fc,
             y = snrna_fc)) + 
  geom_point(aes(color = deg)) + 
  theme_linedraw()
