# Compares this dataset's pseudobulk DESeq2 results (05_tdp43_de_pseudobulk.R)
# with Novartis's TDP-43 knockout Ngn2 neuron DE tables (Ribozero bulk, and
# Flex single-cell and single-nucleus).

suppressMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(ComplexUpset)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

# Matches the cell type tested in the 05_tdp43_de_* scripts.
cell_type <- "L2/3 EN"
cell_type_file <- str_replace_all(cell_type, "[ /]", "_")

results_dir <- "results/sc_results/tdp43neg_vs_novartis_data/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

# The three Novartis tables live in external_data/ (gitignored).
message2("Reading DE tables")

de1 <- read.csv(paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/deseq2.csv"))

# vs Ribozero data --------------------------------------------------------

de2 <- read_xlsx("external_data/TDP43 KO F6 DE Ribozero Ngn2.xlsx",
                 sheet = 1)

str(de1)

de1 <- de1[, c(1, 3, 7)]
de2 <- de2[, c(8, 3, 7)]

colnames(de1) <- c("gene", "nu_lfc", "nu_padj")
colnames(de2) <- c("gene", "novartis_lfc", "novartis_padj")

de <- full_join(de1, de2, by = "gene")

labs <- de %>% 
  na.omit() %>%
  filter(nu_padj < 0.05)

p <- de %>% 
  na.omit() %>%
  mutate(nu_sig = if_else(nu_padj < 0.05, "sig", "ns") %>%
           factor(levels = c("sig", "ns"),
                  labels = c("Yes", "No"))) %>%
  ggplot(aes(x = novartis_lfc,
             y = nu_lfc)) + 
  geom_point(aes(fill = nu_sig),
             shape = 21,
             size = 5) + 
  geom_label_repel(data = labs,
                   aes(label = gene)) +
  # scale_fill_manual(values = c("navy", "orange")) +
  labs(y = "Fold change in TDP-43- neuronal nuclei",
       x = "Fold change in TDP-43 KO Ngn2 cells",
       fill = "Padj < 0.05 in\nTDP-43- neuronal nuclei") +
  theme_linedraw()
ggsave(p,
       filename = paste0(results_dir, "tdp43ko_ribozero_vs_tdp43neg.png"),
       units = "in", dpi = 600,
       height = 8, width = 10)

# ggplot(de2,
#        aes(x = novartis_lfc,
#            y = -log10(novartis_padj))) + 
#   geom_point()

# vs Flex data ------------------------------------------------------------

de3 <- read.csv("external_data/TDP43 KO GW FLEX_DE_results_scRNA copy.csv")
de4 <- read.csv("external_data/TDP43 KO GW FLEX_DE_results_snRNA copy.csv")

de3 <- de3[,c(1, 4, 6)]
de4 <- de4[,c(1, 4, 6)]

colnames(de3) <- c("gene", "novartis_sc_lfc", "novartis_sc_padj")
colnames(de4) <- c("gene", "novartis_sn_lfc", "novartis_sn_padj")

sc_de <- full_join(de1, de3, by = "gene")
sn_de <- full_join(de1, de4, by = "gene")

sc_de <- sc_de %>% 
  mutate(sig = case_when(.default = "ns",
                         nu_padj < 0.05 & abs(nu_lfc) > 0.5 & novartis_sc_padj < 0.05 & abs(novartis_sc_lfc) > 0.5 ~ "both",
                         nu_padj < 0.05 & abs(nu_lfc) > 0.5 & novartis_sc_padj >= 0.05 ~ "nu",
                         nu_padj >= 0.05 & novartis_sc_padj < 0.05 & abs(novartis_sc_lfc) > 0.5 ~ "novartis")) %>% 
  filter(sig != "ns")

p <- sc_de %>%
  na.omit() %>%
  mutate(sig = factor(sig, 
                      levels = c("nu", "novartis", "both"),
                      labels = c("TDP-43- neuronal nuclei", "TDP-43 KO cells", "Both"))) %>%
  # mutate(nu_sig = if_else(nu_padj < 0.05, "sig", "ns") %>%
  #          factor(levels = c("sig", "ns"),
  #                 labels = c("Yes", "No"))) %>%
  ggplot(aes(x = novartis_sc_lfc,
             y = nu_lfc)) + 
  geom_point(aes(fill = sig),
             shape = 21,
             size = 5) + 
  # geom_label_repel(data = labs,
  #                  aes(label = gene)) +
  # scale_fill_manual(values = c("navy", "orange")) +
  labs(y = "Fold change in TDP-43- neuronal nuclei",
       x = "Fold change in TDP-43 KO Ngn2 cells (single-cell)",
       fill = "Padj < 0.05 in:") +
  theme_linedraw()
ggsave(p,
       filename = paste0(results_dir, "tdp43ko_flex_sc_vs_tdp43neg.png"),
       units = "in", dpi = 600,
       height = 8, width = 10)

sn_de <- sn_de %>%
  mutate(sig = case_when(.default = "ns",
                         nu_padj < 0.05 & abs(nu_lfc) > 0.5 & novartis_sn_padj < 0.05 & abs(novartis_sn_lfc) > 0.5 ~ "both",
                         nu_padj < 0.05 & abs(nu_lfc) > 0.5 & novartis_sn_padj >= 0.05 ~ "nu",
                         nu_padj >= 0.05 & novartis_sn_padj < 0.05 & abs(novartis_sn_lfc) > 0.5 ~ "novartis")) %>% 
  filter(sig != "ns")

p <- sn_de %>%
  na.omit() %>%
  mutate(sig = factor(sig, 
                      levels = c("nu", "novartis", "both"),
                      labels = c("TDP-43- neuronal nuclei", "TDP-43 KO nuclei", "Both"))) %>%
  # mutate(nu_sig = if_else(nu_padj < 0.05, "sig", "ns") %>%
  #          factor(levels = c("sig", "ns"),
  #                 labels = c("Yes", "No"))) %>%
  ggplot(aes(x = novartis_sn_lfc,
             y = nu_lfc)) + 
  geom_point(aes(fill = sig),
             shape = 21,
             size = 5) + 
  # geom_label_repel(data = labs,
  #                  aes(label = gene)) +
  # scale_fill_manual(values = c("navy", "orange")) +
  labs(y = "Fold change in TDP-43- neuronal nuclei",
       x = "Fold change in TDP-43 KO Ngn2 cells (single-nuclei)",
       fill = "Padj < 0.05 in:") +
  theme_linedraw()
ggsave(p,
       filename = paste0(results_dir, "tdp43ko_flex_sn_vs_tdp43neg.png"),
       units = "in", dpi = 600,
       height = 8, width = 10)

# upset plot --------------------------------------------------------------

list <- list()
list[[1]] <- de1 %>%
  filter(abs(nu_lfc) > 0.5 & nu_padj < 0.05) %>% 
  pull(gene)
list[[2]] <- de2 %>% 
  filter(abs(novartis_lfc) > 0.5 & novartis_padj < 0.05) %>% 
  pull(gene)
list[[3]] <- de3 %>% 
  filter(abs(novartis_sc_lfc) > 0.5 & novartis_sc_padj < 0.05) %>% 
  pull(gene)
list[[4]] <- de4 %>%
  filter(abs(novartis_sn_lfc) > 0.5 & novartis_sn_padj < 0.05) %>% 
  pull(gene)

names(list) <- c("TDP-43- nuclei", "TDP-43 KO Ribozero", "TDP-43 KO sc Flex", "TDP-43 KO sn Flex")

genes <- unique(list_c(list))

intersect_df <- data.frame(gene = genes)
for (i in seq_along(list)){
  
  col <- names(list)[i]
  
  intersect_df <- intersect_df %>%
    mutate(x = if_else(gene %in% list[[i]], T, F)) %>%
    dplyr::rename(!!sym(col) := "x")
}

p_upset <- upset(intersect_df,
      names(list),
      height_ratio = 0.9,
      # width_ratio = 1.5,
      encode_sets = F,
      # base_annotations=list(
      #   'Intersection size'=intersection_size(counts=FALSE)
      # ),
      themes = upset_modify_themes(
        list(
          'Intersection size' = theme(axis.text = element_text(color = "black"),
                                      axis.ticks.y = element_line(),
                                      panel.border = element_rect(color = "black", fill = "transparent")),
          'intersections_matrix' = theme(axis.text = element_text(color = "black"),
                                         axis.title = element_blank(),
                                         panel.border = element_rect(color = "black", fill = "transparent")),
          'overall_sizes' = theme(panel.border = element_rect(color = "black", fill = "transparent"))
        )
      ),
      set_sizes = upset_set_size() + 
        scale_y_continuous(labels = function(x) x / 1000) + 
        labs(y = "Set size (1000s)")
)

sc_de_labs <- sc_de %>%
  filter(sig == "both") %>%
  mutate(fc1 = abs(nu_lfc),
         fc2 = abs(novartis_sc_lfc)) %>% 
  mutate(rank = fc1 * fc2) %>% 
  arrange(desc(rank)) %>%
  slice_head(n = 30)

p_sc_de <- sc_de %>%
  na.omit() %>%
  mutate(sig = factor(sig, 
                      levels = c("nu", "novartis", "both"),
                      labels = c("TDP-43- neuronal nuclei", "TDP-43 KO cells", "Both"))) %>%
  # mutate(nu_sig = if_else(nu_padj < 0.05, "sig", "ns") %>%
  #          factor(levels = c("sig", "ns"),
  #                 labels = c("Yes", "No"))) %>%
  ggplot(aes(x = novartis_sc_lfc,
             y = nu_lfc)) + 
  geom_point(aes(fill = sig),
             shape = 21,
             size = 5) + 
  geom_label_repel(data = sc_de_labs,
                   aes(label = gene),
                   size = 3) +
  ggtitle("TDP-43- nuclei vs TDP-43 KO cells") +
  labs(y = "Fold change in TDP-43- neuronal nuclei",
       x = "Fold change in TDP-43 KO Ngn2 cells (single-cell)",
       fill = "Padj < 0.05 in:") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5))

sn_de_labs <- sn_de %>%
  filter(sig == "both") %>%
  mutate(fc1 = abs(nu_lfc),
         fc2 = abs(novartis_sn_lfc)) %>% 
  mutate(rank = fc1 * fc2) %>% 
  arrange(desc(rank)) %>%
  slice_head(n = 30)

p_sn_de <- sn_de %>%
  na.omit() %>%
  mutate(sig = factor(sig, 
                      levels = c("nu", "novartis", "both"),
                      labels = c("TDP-43- neuronal nuclei", "TDP-43 KO nuclei", "Both"))) %>%
  # mutate(nu_sig = if_else(nu_padj < 0.05, "sig", "ns") %>%
  #          factor(levels = c("sig", "ns"),
  #                 labels = c("Yes", "No"))) %>%
  ggplot(aes(x = novartis_sn_lfc,
             y = nu_lfc)) + 
  geom_point(aes(fill = sig),
             shape = 21,
             size = 5) + 
  geom_label_repel(data = sn_de_labs,
                   aes(label = gene),
                   size = 3) +
  # scale_fill_manual(values = c("navy", "orange")) +
  ggtitle("TDP-43- nuclei vs TDP-43 KO nuclei") +
  labs(y = "Fold change in TDP-43- neuronal nuclei",
       x = "Fold change in TDP-43 KO Ngn2 cells (single-nucleus)",
       fill = "Padj < 0.05 in:") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5))

p <- free(p_upset) / (p_sc_de | p_sn_de) 
ggsave(p,
       filename = paste0(results_dir, "tdp43neg_vs_tdp43ko_flex.png"),
       units = "in", dpi = 600,
       height = 10, width = 15)
