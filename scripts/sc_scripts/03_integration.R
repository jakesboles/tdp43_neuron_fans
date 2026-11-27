# Normalizes and integrates every QC-filtered sample (Harmony, across
# samples), clusters, and finds cluster markers. Rebuilds the object from
# 01_qc.R's BPCells counts, adds 02_doubletfinder.R's calls, and saves
# metadata + reductions for 04_annotation.R.

suppressMessages({
  library(tidyverse)
  library(scCustomize)
  library(Seurat)
  library(BPCells)
  library(reticulate)
  library(dittoSeq)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

results_dir <- "results/sc_results/03_integration/"
data_dir <- "data/sc_data/03_integration/"

for (i in c(results_dir, paste0(results_dir, "pca_loadings/"), data_dir)){
  dir.create(i,
             showWarnings = F,
             recursive = T)
}

# Load QC-filtered counts ---------------------------------------------------
message2("Loading QC-filtered counts")

meta <- readRDS("data/sc_data/01_qc/metadata.rds")
counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")

# JCO_Four() only has 4 colors, so fall back to a larger palette once there
# are more than 4 samples (i.e. more than one pool).
n_samples <- length(unique(obj$orig.ident))
if (n_samples <= 4){
  sample_colors <- JCO_Four()
} else {
  sample_colors <- DiscretePalette_scCustomize(num_colors = n_samples,
                                               palette = "polychrome")
}

# Add DoubletFinder info to samples ---------------------------------------
message2("Adding DoubletFinder calls")

doubletfinder_files <- list.files("results/sc_results/02_doubletfinder",
                                  pattern = "_meta\\.csv$",
                                  recursive = F,
                                  full.names = T)

# Fails loudly if any 02_doubletfinder.R array task didn't finish, instead
# of silently leaving those samples' DF.unadj as NA.
df_samples <- str_remove(basename(doubletfinder_files), "_meta\\.csv$")
missing <- setdiff(unique(obj$orig.ident), df_samples)
if (length(missing) > 0){
  stop(paste0(length(missing), " sample(s) missing from ",
              "results/sc_results/02_doubletfinder/ -- 02_doubletfinder.R may ",
              "not have finished (or failed) for: ",
              paste(missing, collapse = ", ")))
}

doubletfinder_tab <- map(doubletfinder_files,
                         read.csv) %>%
  list_rbind()

table(doubletfinder_tab$DF.unadj)

doubletfinder_tab <- doubletfinder_tab %>%
  column_to_rownames(var = "X") %>%
  dplyr::select(c("DF.unadj"))

obj <- AddMetaData(obj,
                   doubletfinder_tab)

obj <- Split_Layers(obj,
                    split.by = "orig.ident")

# Process data ------------------------------------------------------------
message2("Normalizing, scaling, and running PCA")

obj <- NormalizeData(obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

p <- ElbowPlot(obj)
ggsave(p,
       filename = paste0(results_dir, "pca_elbow.png"),
       units = "in", dpi = 600,
       height = 5, width = 5,
       bg = "white")

# One PNG per PC, replacing Iterate_PC_Loading_Plots(), which only writes a
# single multi-page PDF.
for (i in 1:50){
  p <- PC_Plotting(obj,
                   dim_number = i)
  ggsave(p,
         filename = paste0(results_dir, "pca_loadings/PC_", i, ".png"),
         units = "in", dpi = 300,
         height = 8, width = 6,
         bg = "white")
}

message2("Computing unintegrated UMAP")

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

p <- DimPlot_scCustom(obj,
                      group.by = "orig.ident",
                      reduction = "raw_umap",
                      colors_use = sample_colors)
ggsave(p,
       filename = paste0(results_dir, "raw_umap_samples.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

message2("Integrating samples using Harmony")

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

p <- DimPlot_scCustom(obj,
                      group.by = "orig.ident",
                      reduction = "integrated_umap",
                      colors_use = sample_colors)
ggsave(p,
       filename = paste0(results_dir, "integrated_umap_samples.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      group.by = "tissue",
                      reduction = "integrated_umap") +
  DimPlot_scCustom(obj,
                   group.by = "tdp",
                   reduction = "integrated_umap",
                   colors_use = c("magenta3", "springgreen2"))
ggsave(p,
       filename = paste0(results_dir, "integrated_umap_split.png"),
       units = "in", dpi = 600,
       height = 5, width = 12)

message2("Clustering")

obj <- obj %>%
  FindClusters(algorithm = 4,
               method = "igraph",
               cluster.name = "leiden_cluster",
               resolution = 1)

p <- DimPlot_scCustom(obj,
                      reduction = "integrated_umap")
ggsave(p,
       filename = paste0(results_dir, "integrated_umap_clusters.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

# dittoSeq sorts cluster ids as text (1, 10, 11, ..., 2, ...); this puts
# them back in numeric order for however many clusters there are, instead
# of the hard-coded c(1, 11:18, 2:10), which only fit 18 clusters.
cluster_levels <- levels(factor(as.character(obj$leiden_cluster)))
cluster_order <- order(as.numeric(cluster_levels))

p <- dittoBarPlot(obj,
                  group.by = "leiden_cluster",
                  var = "orig.ident",
                  color.panel = sample_colors,
                  x.reorder = cluster_order) +
  scale_y_continuous(expand = c(0, 0))
ggsave(p,
       filename = paste0(results_dir, "cluster_by_sample_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 8)

# Renamed: this previously saved to cluster_by_sample_bars.png too,
# overwriting the plot above.
# NavyAndOrange() only has 2 colors, so fall back to dittoSeq's default
# palette if new pools bring more than 2 tissue values.
if (length(unique(obj$tissue)) <= 2){
  tissue_colors <- NavyAndOrange()
} else {
  tissue_colors <- dittoColors()
}

p <- dittoBarPlot(obj,
                  group.by = "leiden_cluster",
                  var = "tissue",
                  color.panel = tissue_colors,
                  x.reorder = cluster_order) +
  scale_y_continuous(expand = c(0, 0)) +
  dittoBarPlot(obj,
               group.by = "leiden_cluster",
               var = "tdp",
               color.panel = c("magenta3", "springgreen2"),
               x.reorder = cluster_order) +
  scale_y_continuous(expand = c(0, 0))
ggsave(p,
       filename = paste0(results_dir, "cluster_by_tissue_tdp_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 15)

obj <- JoinLayers(obj)

# Find markers --------------------------------------------------------------
message2("Finding cluster markers")

# FindAllMarkers() tests gene by gene, but BPCells data here is stored
# cell-major; transposing once up front (to a temp dir) avoids Seurat
# re-transposing small chunks throughout the test loop.
obj[["RNA"]]$data <- transpose_storage_order(obj[["RNA"]]$data)

Idents(obj) <- "leiden_cluster"

markers <- FindAllMarkers(obj)

write.csv(markers,
          file = paste0(results_dir, "cluster_markers.csv"))

# Save ----------------------------------------------------------------------
message2("Saving metadata and reductions")

# No matrix is written here: 04_annotation.R rebuilds counts from 01_qc.R's
# BPCells matrix and re-normalizes (cheap and deterministic).
saveRDS(obj@meta.data,
        file = paste0(data_dir, "metadata.rds"))

saveRDS(obj[["integrated"]],
        file = paste0(data_dir, "integrated.rds"))

saveRDS(obj[["integrated_umap"]],
        file = paste0(data_dir, "integrated_umap.rds"))
