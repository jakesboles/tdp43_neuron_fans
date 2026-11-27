library(tidyverse)
library(scCustomize)
library(Seurat)
library(reticulate)
library(dittoSeq)

setwd("/projects/b1169/boles/JSB149")

plots_dir <- "plots/03_integration/"
dir.create(plots_dir,
           recursive = T)

data_dir <- "data/03_integration/"
dir.create(data_dir,
           recursive = T)

tab_dir <- "tab_data/03_integration/"
dir.create(tab_dir,
           recursive = T)

obj <- readRDS("data/01_qc/obj.rds")

obj <- Split_Layers(obj,
                    split.by = "orig.ident")

# Add DoubletFinder info to samples ---------------------------------------

doubletfinder_files <- list.files("tab_data/02_doubletfinder",
                                  recursive = F,
                                  full.names = T)

doubletfinder_tabs <- map(doubletfinder_files,
                          read.csv)

doubletfinder_tab <- list_rbind(doubletfinder_tabs)

table(doubletfinder_tab$DF.unadj)

doubletfinder_tab <- doubletfinder_tab %>%
  column_to_rownames(var = "X") %>%
  dplyr::select(c("DF.unadj"))

obj <- AddMetaData(obj,
                   doubletfinder_tab)

# Process data ------------------------------------------------------------

obj <- NormalizeData(obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(obj)
ggsave(filename = paste0(plots_dir, "pca_elbow.png"),
       units = "in", dpi = 600,
       height = 5, width = 5,
       bg = "white")

Iterate_PC_Loading_Plots(obj,
                         dims_plot = 50,
                         file_path = plots_dir,
                         file_name = "pca_loadings")

obj <- obj %>%
  FindNeighbors(reduction = "pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)

obj <- obj %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "raw_umap",
          return.model = F)

DimPlot_scCustom(obj,
                 group.by = "orig.ident",
                 reduction = "raw_umap",
                 colors_use = JCO_Four())
ggsave(filename = paste0(plots_dir, "raw_umap_samples.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       assay = "RNA",
                       npcs = 10,
                       k.anchor = 20,
                       new.reduction = "integrated")

obj <- obj %>%
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
                compute.SNN = T)

obj <- obj %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "integrated_umap",
          return.model = F)

DimPlot_scCustom(obj,
                 group.by = "orig.ident",
                 reduction = "integrated_umap",
                 colors_use = JCO_Four())
ggsave(filename = paste0(plots_dir, "integrated_umap_samples.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

DimPlot_scCustom(obj,
                 group.by = "tissue",
                 reduction = "integrated_umap") + 
DimPlot_scCustom(obj,
                 group.by = "tdp",
                 reduction = "integrated_umap",
                 colors_use = c("magenta3", "springgreen2"))
ggsave(filename = paste0(plots_dir, "integrated_umap_split.png"),
       units = "in", dpi = 600,
       height = 5, width = 12)

obj <- obj %>%
  FindClusters(algorithm = 4,
               method = "igraph",
               cluster.name = "leiden_cluster",
               resolution = 1)

DimPlot_scCustom(obj,
                 reduction = "integrated_umap")
ggsave(filename = paste0(plots_dir, "integrated_umap_clusters.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

dittoBarPlot(obj,
             group.by = "leiden_cluster",
             var = "orig.ident",
             color.panel = JCO_Four(),
             x.reorder = c(1, 11:18, 2:10)) + 
  scale_y_continuous(expand = c(0, 0))
ggsave(filename = paste0(plots_dir, "cluster_by_sample_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 8)

dittoBarPlot(obj,
             group.by = "leiden_cluster",
             var = "tissue",
             color.panel = NavyAndOrange(),
             x.reorder = c(1, 11:18, 2:10)) + 
  scale_y_continuous(expand = c(0, 0))+
  dittoBarPlot(obj,
               group.by = "leiden_cluster",
               var = "tdp",
               color.panel = c("magenta3", "springgreen2"),
               x.reorder = c(1, 11:18, 2:10)) + 
  scale_y_continuous(expand = c(0, 0))
ggsave(filename = paste0(plots_dir, "cluster_by_sample_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 15)

obj <- JoinLayers(obj)

markers <- FindAllMarkers(obj)

write.csv(markers,
          file = paste0(tab_dir, "cluster_markers.csv"))

saveRDS(obj,
        file = paste0(data_dir, "obj.rds"))