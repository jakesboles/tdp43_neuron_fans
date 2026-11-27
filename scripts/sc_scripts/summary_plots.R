# Summary figure: cell types, nuclear TDP-43 status, TDP-43 status by cell
# type, and the pseudobulk DESeq2 volcano for one cell type.

suppressMessages({
  library(Seurat)
  library(BPCells)
  library(scCustomize)
  library(tidyverse)
  library(dittoSeq)
  library(ggrepel)
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

results_dir <- "results/sc_results/summary_plots/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

# dittoSeq orders groups as text; this returns the x.reorder indices that
# put whichever of `wanted` are present in that order (see 04_annotation.R).
ditto_order <- function(x, wanted){
  lv <- levels(factor(as.character(x)))
  match(wanted[wanted %in% lv], lv)
}

message2("Loading metadata, UMAP, and DE results")

meta <- readRDS("data/sc_data/04_annotation/metadata.rds")
counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
obj[["integrated_umap"]] <- readRDS("data/sc_data/04_annotation/integrated_umap.rds")

de <- read.csv(paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/deseq2.csv"))

colnames(obj@meta.data)

p1 <- DimPlot_scCustom(obj,
                 group.by = "label2",
                 reduction = "integrated_umap") + 
  ggtitle("Neuronal types")

p2 <- DimPlot_scCustom(obj,
                 group.by = "tdp",
                 reduction = "integrated_umap") + 
  ggtitle("Nuclear TDP-43")

p3 <- dittoBarPlot(obj,
             var = "tdp",
             group.by = "label2",
             x.reorder = ditto_order(obj$label2, levels(obj$label2)),
             color.panel = NavyAndOrange()) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = rev) +
  ggtitle("Nuclear TDP-43 by cell type") + 
  theme_linedraw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) + 
  coord_flip()

de <- de %>% 
  mutate(sig = case_when(padj < 0.05 & LFC > 0.25 ~ "up",
                       padj < 0.05 & LFC < -0.25  ~ "down",
                       .default = "ns"),
         rank = LFC * -log10(padj))

labs1 <- de %>%
  filter(padj < 0.05 & 
           str_detect(X, "ENSG|LINC", negate = T)) %>% 
  arrange(desc(rank)) %>% 
  slice_head(n = 20)
labs2 <- de %>% 
  filter(padj < 0.05 & 
           str_detect(X, "ENSG|LINC", negate = T)) %>% 
  arrange(rank) %>% 
  slice_head(n = 20)

labs <- rbind(labs1, labs2)

p4 <- de %>% 
  ggplot(aes(x = LFC,
             y = -log10(padj),
             color = sig)) + 
  geom_point(show.legend = F) + 
  geom_label_repel(data = labs,
                  aes(label = X),
                  show.legend = F) +
  scale_color_manual(values = c("orange", "gray60", "navy")) +
  labs(x = "log2FoldChange") +
  ggtitle(paste0("DE in TDP-43- ", cell_type, "s")) +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


design <- "
AABB
CCCC
DDDD
"

p <- p1 + p2 + p3 + p4 + 
  plot_layout(design = design,
              heights = c(1, 0.5, 1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 32, face = "bold"))
ggsave(p,
       filename = paste0(results_dir, "summary.png"),
       units = "in", dpi = 600,
       height = 10, width = 12)
