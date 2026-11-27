library(tidyverse)
library(scCustomize)
library(Seurat)
library(dittoSeq)
library(patchwork)

setwd("/projects/b1169/boles/JSB149")

plots_dir <- "plots/04_annotation/"
tab_dir <- "tab_data/04_annotation/"
data_dir <- "data/04_annotation/"

for (i in c(plots_dir, tab_dir, data_dir)){
  dir.create(i,
             showWarnings = F,
             recursive = T)
}

obj <- readRDS("data/03_integration/obj.rds")

markers <- read.csv("tab_data/03_integration/cluster_markers.csv")

top <- markers %>% 
  filter(pct.1 > 0.3 & 
           str_detect(gene, "ENSG|LINC", negate = T)) %>% 
  Extract_Top_Markers(num_genes = 3,
                      make_unique = T,
                      named_vector = F)

top

dittoDotPlot(obj,
             vars = top,
             group.by = "leiden_cluster")

markers %>%
  filter(cluster == 14 & 
           pct.1 > 0.3) %>% 
  arrange(desc(avg_log2FC)) %>%
  head(20)

FeaturePlot_scCustom(obj,
                     features = c("SLC17A6", "GAD1"),
                     reduction = "integrated_umap")

DimPlot_scCustom(obj,
                 group.by = "DF.unadj",
                 reduction = "integrated_umap")

markers %>% 
  filter(gene == "CLDN5") %>%
  arrange(desc(avg_log2FC))

obj@meta.data <- obj@meta.data %>% 
  mutate(label1 = case_when(leiden_cluster %in% 14 ~ "Contamination",
                            leiden_cluster %in% c(6, 9, 11, 12) ~ "Inhibitory neuron",
                            leiden_cluster %in% c(1, 2, 3, 4, 5, 7, 8, 10, 13, 15, 16, 17, 18) ~ "Excitatory neuron"))

DimPlot_scCustom(obj,
                 group.by = "label1",
                 reduction = "integrated_umap",
                 colors_use = c("cyan3", "yellow3", "magenta3"))
ggsave(filename = paste0(plots_dir, "broad_labels.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

# Inhibitory neurons ------------------------------------------------------

inh <- subset(obj,
             subset = label1 == "Inhibitory neuron") %>% 
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "integrated_umap",
          return.model = F) %>%
  FindClusters(algorithm = 4,
               method = "igraph",
               cluster.name = "leiden_cluster",
               resolution = 1)

DimPlot_scCustom(inh,
                 reduction = "integrated_umap",
                 pt.size = 1.5)
ggsave(filename = paste0(plots_dir, "inhibitory_neuron_clusters.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

inh_markers <- FindAllMarkers(inh)
write.csv(inh_markers,
          file = paste0(tab_dir, "inhibitory_neuron_markers.csv"))

inh_top <- inh_markers %>% 
  filter(pct.1 > 0.3 & 
           str_detect(gene, "ENSG|LINC", negate = T)) %>% 
  Extract_Top_Markers(5,
                      make_unique = T,
                      named_vector = F)

dittoDotPlot(inh,
             vars = inh_top,
             group.by = "leiden_cluster")
ggsave(filename = paste0(plots_dir, "inh_markers.png"),
       units = "in", dpi = 600,
       height = 5, width = 10)

{gene <- "VIP"
  FeaturePlot_scCustom(inh,
                       features = gene,
                       reduction = "integrated_umap") + 
    VlnPlot_scCustom(inh,
                     features = gene)
}

# query <- FindMarkers(inh,
#                      ident.1 = c(3, 6))
# 
# query %>%
#   filter(pct.1 > 0.3) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   head(20)

inh_markers %>% 
  filter(cluster == 10 & 
           pct.1 > 0.3 & 
           str_detect(gene, "LINC|ENSG", negate = T)) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(20)

dittoBarPlot(inh,
             var = "DF.unadj",
             group.by = "leiden_cluster")

inh@meta.data <- inh@meta.data %>% 
  mutate(label2 = case_when(leiden_cluster %in% c(4, 8, 9) ~ "LAMP5 IN",
                            leiden_cluster %in% c(1, 5, 10) ~ "VIP IN",
                            leiden_cluster %in% c(2, 11) ~ "CXCL14 IN",
                            leiden_cluster %in% c(3, 12) ~ "PVALB IN",
                            leiden_cluster %in% 7 ~ "SST IN",
                            leiden_cluster == 6 ~ "Doublet"))


DimPlot_scCustom(inh,
                 group.by = "label2",
                 reduction = "integrated_umap",
                 colors_use = Dark2_Pal())
ggsave(filename = paste0(plots_dir, "inhibitory_neuron_labels.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Excitatory neurons ------------------------------------------------------


ex <- subset(obj,
             subset = label1 == "Excitatory neuron") %>% 
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "integrated_umap",
          return.model = F) %>%
  FindClusters(algorithm = 4,
               method = "igraph",
               cluster.name = "leiden_cluster",
               resolution = 1.5)

DimPlot_scCustom(ex,
                 reduction = "integrated_umap",
                 pt.size = 1.5)
ggsave(filename = paste0(plots_dir, "excitatory_neuron_clusters.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

dittoBarPlot(ex,
             var = "DF.unadj",
             group.by = "leiden_cluster")

dittoBarPlot(ex,
             var = "tdp",
             group.by = "leiden_cluster")

ex_markers <- FindAllMarkers(ex)
write.csv(ex_markers,
          file = paste0(tab_dir, "excitatory_neuron_markers.csv"))

ex_top <- ex_markers %>% 
  filter(pct.1 > 0.3 & 
           str_detect(gene, "ENSG|LINC", negate = T)) %>% 
  Extract_Top_Markers(5,
                      make_unique = T,
                      named_vector = F)

dittoDotPlot(ex,
             vars = ex_top,
             group.by = "leiden_cluster")
ggsave(filename = paste0(plots_dir, "ex_markers.png"),
       units = "in", dpi = 600,
       height = 5, width = 10)

{gene <- "HTR2C"
  FeaturePlot_scCustom(ex,
                       features = gene,
                       reduction = "integrated_umap") + 
    VlnPlot_scCustom(ex,
                     features = gene)
}

ex_markers %>%
  filter(cluster == 14 & 
           pct.1 > 0.3 & 
           str_detect(gene, "LINC|ENSG", negate = T)) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(20)

ex@meta.data <- ex@meta.data %>%
  mutate(label2 = case_when(leiden_cluster %in% c(1, 2, 3, 5, 7, 8, 9, 10, 11, 14, 15) ~ "L2/3 EN",
                            leiden_cluster %in% c(6, 17) ~ "L4 EN",
                            leiden_cluster %in% c(12) ~ "L4/5 EN",
                            leiden_cluster %in% c(18, 16, 4, 13, 19) ~ "L5/6 EN"))

DimPlot_scCustom(ex,
                 reduction = "integrated_umap",
                 group.by = "label2",
                 colors_use = Dark2_Pal())

ex_meta <- ex@meta.data
inh_meta <- inh@meta.data
full_meta <- obj@meta.data
new_meta <- rbind(ex_meta,
                  inh_meta)
obj <- AddMetaData(obj,
                   new_meta)
DimPlot_scCustom(obj,
                 group.by = "label2",
                 reduction = "integrated_umap")
final <- subset(obj,
                subset = label2 != "Contamination" & label2 != "Doublet")

final <- Split_Layers(final,
                      split.by = "orig.ident")

final <- final %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  IntegrateLayers(method = "HarmonyIntegration",
                  orig.reduction = "pca",
                  assay = "RNA",
                  npcs = 10,
                  k.anchor = 20,
                  new.reduction = "integrated") %>%
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "integrated",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "integrated_umap",
          return.model = F)

final$label2 <- factor(final$label2,
                       levels = c("L2/3 EN", "L4 EN", "L4/5 EN", "L5/6 EN", "CXCL14 IN", "LAMP5 IN", "PVALB IN", "SST IN", "VIP IN"))

p1 <- DimPlot_scCustom(final,
                       group.by = "label2",
                       reduction = "integrated_umap",
                       colors_use = "glasbey")
p2 <- dittoBarPlot(final,
                   var = "tdp",
                   group.by = "label2",
                   main = "Nuclear TDP content by cell type",
                   color.panel = c("navy", "orange"),
                   x.reorder = c(2, 3, 4, 5, 1, 6, 7, 8, 9)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

p3 <- DimPlot_scCustom(final,
                       group.by = "tdp",
                       reduction = "integrated_umap",
                       pt.size = 1) + 
  theme(plot.title = element_blank())

p1 + p3 + p2 + 
  plot_layout(ncol = 3)

# save barcodes for Jackie
write.table(colnames(final),
          file = paste0(tab_dir, "filtered_barcodes.csv"),
          row.names = F,
          quote = F,
          col.names = F)

# final@meta.data %>%
#   dplyr::select(orig.ident) %>% 
#   rownames_to_column(var = "barcode") %>% 
#   mutate(barcode = str_split_i(barcode, "_", i = 3)) %>%
#   dplyr::rename("alias" = "orig.ident") %>%
#   write.csv()

saveRDS(final,
        file = paste0(data_dir, "obj.rds"))
