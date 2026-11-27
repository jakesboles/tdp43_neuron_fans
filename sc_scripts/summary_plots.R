library(Seurat)
library(scCustomize)
library(tidyverse)
library(dittoSeq)
library(ggrepel)
library(patchwork)

setwd("/projects/b1169/boles/JSB149")

obj <- readRDS("data/04_annotation/obj.rds")

de <- read.csv("tab_data/05_tdp43_de/deseq2.csv")

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
             x.reorder = c(2, 3, 4, 1, 5, 6, 7, 8),
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
  ggtitle("DE in TDP-43- L2/3 neurons") +
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
       filename = "plots/jsb149_summary.png",
       units = "in", dpi = 600,
       height = 10, width = 12)
