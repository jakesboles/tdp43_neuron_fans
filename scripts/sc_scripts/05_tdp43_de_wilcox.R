# Wilcoxon (FindMarkers) DE, TDPneg vs TDPpos, within one cell type, plus
# the lab-meeting overview figure. Rebuilds counts from 01_qc.R's BPCells
# matrix using 04_annotation.R's labels and UMAP.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(scCustomize)
  library(patchwork)
  library(ggrepel)
  library(dittoSeq)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

# NEEDS INPUT: only L2/3 ENs are tested, as before (see
# 05_tdp43_de_mast.R).
cell_type <- "L2/3 EN"
cell_type_file <- str_replace_all(cell_type, "[ /]", "_")

results_dir <- paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/")
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

# dittoSeq orders groups as text; this returns the x.reorder indices that
# put whichever of `wanted` are present in that order (see 04_annotation.R).
ditto_order <- function(x, wanted){
  lv <- levels(factor(as.character(x)))
  match(wanted[wanted %in% lv], lv)
}

# Load 04_annotation.R's object ---------------------------------------------
message2("Loading counts, metadata, and UMAP")

meta <- readRDS("data/sc_data/04_annotation/metadata.rds")
counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

obj0 <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
obj0[["integrated_umap"]] <- readRDS("data/sc_data/04_annotation/integrated_umap.rds")

# One cell type fits in memory as a plain matrix, which avoids BPCells'
# storage-order slowdown in FindMarkers().
obj <- subset(obj0, label2 == cell_type)
obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, "dgCMatrix")
obj <- NormalizeData(obj)

Idents(obj) <- "tdp"

message2("Running Wilcoxon DE")

de <- FindMarkers(obj,
                  ident.1 = "TDPneg",
                  ident.2 = "TDPpos",
                  min.pct = 0.1)

de %>%
  arrange(p_val_adj) %>%
  rownames_to_column(var = "gene") %>%
  filter(str_detect(gene, "ENSG|LINC", negate = T)) %>%
  head(30)

de["ATXN2",]

p <- VlnPlot_scCustom(obj,
                      features = "ATXN2")
ggsave(p,
       filename = paste0(results_dir, "wilcox_ATXN2_vln.png"),
       units = "in", dpi = 600,
       height = 4, width = 4)

labs1 <- de %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter( p_val_adj < 0.05 &
           str_detect(gene, "ENSG|LINC", negate = T)) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30)

labs2 <- de %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 0.05 &
            str_detect(gene, "ENSG|LINC", negate = T)) %>%
  arrange(avg_log2FC) %>%
  slice_head(n =30)

labs3 <- de %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% c("UNC13A", "STMN2", "HDGFL2", "KCNG2", "PKBP"))

labs <- rbind(labs1, labs2) %>%
  rbind(labs3)

p <- de %>%
  mutate(sig = case_when(avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
                         avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
                         .default = "ns")) %>%
  mutate(sig = factor(sig,
                      levels = c("Up", "Down", "ns"))) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj))) +
  geom_point(aes(color = sig)) +
  geom_label_repel(data = labs,
                   aes(label = gene)) +
  scale_color_manual(values = c("navy", "orange", "gray60")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(p,
       filename = paste0(results_dir, "wilcox_volcano.png"),
       units = "in", dpi = 600,
       height = 7, width = 8)

de %>%
  write.csv(paste0(results_dir, "wilcox.csv"))

# Extra figs for lab meeting ----------------------------------------------
message2("Making lab meeting figure")

p1 <- DimPlot_scCustom(obj0,
                 group.by = "label2",
                 reduction = "integrated_umap",
                 colors_use = "glasbey") +
  NoAxes() +
  theme(plot.title = element_blank(),
        legend.text = element_text(size = 16))

# label2 is saved as a factor in 04_annotation.R's intended order; this
# replaces the hard-coded c(2, 3, 4, 1, 5, 6, 7, 8), which no longer
# matched the 9 labels 04 assigns.
p2 <- dittoBarPlot(obj0,
                   group.by = "label2",
                   var = "tdp",
             color.panel = NavyAndOrange(),
             main = "TDP-43 mislocalization",
             x.reorder = ditto_order(obj0$label2, levels(obj0$label2))) +
  scale_y_continuous(expand = c(0, 0)) +
  # theme_bw(base_size = 16) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

p3 <- de %>%
  mutate(sig = case_when(avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
                         avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
                         .default = "ns")) %>%
  mutate(sig = factor(sig,
                      levels = c("Up", "Down", "ns"))) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj))) +
  geom_point(aes(color = sig)) +
  # geom_label_repel(data = labs,
                   # aes(label = gene)) +
  scale_color_manual(values = c("navy", "orange", "gray60")) +
  ggtitle(paste0("TDP-43- vs TDP-43+\n", cell_type, "s")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

p <- (free(p1) + p2) / p3
ggsave(p,
       filename = paste0(results_dir, "wilcox_lab_meeting.png"),
       units = "in", dpi = 600,
       height = 10, width = 12,
       bg = "white")
