library(tidyverse)
library(Seurat)
library(scCustomize)
library(patchwork)
library(ggrepel)
library(dittoSeq)

setwd("/projects/b1169/boles/JSB149")

dir.create("tab_data/05_tdp43_de",
           showWarnings = F,
           recursive = T)

obj0 <- readRDS("data/04_annotation/obj.rds")

obj <- obj0 %>%
  subset(label2 == "L2/3 EN")

Idents(obj) <- "tdp"

obj <- JoinLayers(obj)

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

VlnPlot_scCustom(obj,
                 features = "ATXN2")

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

de %>%
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

de %>%
  write.csv("tab_data/05_tdp43_de/wilcox.csv")

# Extra figs for lab meeting ----------------------------------------------

p1 <- DimPlot_scCustom(obj0,
                 group.by = "label2",
                 reduction = "integrated_umap",
                 colors_use = "glasbey") + 
  NoAxes() +
  theme(plot.title = element_blank(),
        legend.text = element_text(size = 16))

p2 <- dittoBarPlot(obj0,
                   group.by = "label2",
                   var = "tdp",
             color.panel = NavyAndOrange(),
             main = "TDP-43 mislocalization",
             x.reorder = c(2, 3, 4, 1, 5, 6, 7, 8)) + 
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
  ggtitle("TDP-43- vs TDP-43+\nL2/3 ENs") +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

(free(p1) + p2) / p3
