# Interactive annotation: labels 03_integration.R's clusters as excitatory
# or inhibitory neurons, subclusters each to assign subtypes (label2), drops
# contamination and doublets, then re-integrates the cleaned object. Rebuilds
# the object from 01_qc.R's BPCells counts plus 03_integration.R's metadata
# and reductions.
#
# NEEDS INPUT: every cluster -> label mapping below (case_when blocks) was
# made by hand for pool 149's clustering. Cluster numbers change whenever
# 03_integration.R is rerun on different data (e.g. once more pools are
# added), so these maps must be redone from the new marker plots before the
# labels mean anything.

suppressMessages({
  library(tidyverse)
  library(scCustomize)
  library(Seurat)
  library(BPCells)
  library(dittoSeq)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

results_dir <- "results/sc_results/04_annotation/"
data_dir <- "data/sc_data/04_annotation/"

for (i in c(results_dir, data_dir)){
  dir.create(i,
             showWarnings = F,
             recursive = T)
}

# dittoSeq orders groups as text; this returns the x.reorder indices that
# put whichever of `wanted` are present in that order, so the bar order
# doesn't depend on a hard-coded index vector of a fixed length.
ditto_order <- function(x, wanted){
  lv <- levels(factor(as.character(x)))
  match(wanted[wanted %in% lv], lv)
}

# Load 03_integration.R's object --------------------------------------------
message2("Loading counts, metadata, and reductions")

meta <- readRDS("data/sc_data/03_integration/metadata.rds")
counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
obj <- NormalizeData(obj)
obj[["integrated"]] <- readRDS("data/sc_data/03_integration/integrated.rds")
obj[["integrated_umap"]] <- readRDS("data/sc_data/03_integration/integrated_umap.rds")

Idents(obj) <- "leiden_cluster"

markers <- read.csv("results/sc_results/03_integration/cluster_markers.csv")

# Broad labels --------------------------------------------------------------
message2("Plotting broad markers")

top <- markers %>%
  filter(pct.1 > 0.3 &
           str_detect(gene, "ENSG|LINC", negate = T)) %>%
  Extract_Top_Markers(num_genes = 3,
                      make_unique = T,
                      named_vector = F)

top

p <- dittoDotPlot(obj,
                  vars = top,
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "cluster_top_markers.png"),
       units = "in", dpi = 600,
       height = 5, width = 14,
       bg = "white")

markers %>%
  filter(cluster == 14 &
           pct.1 > 0.3) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

p <- FeaturePlot_scCustom(obj,
                          features = c("SLC17A6", "GAD1"),
                          reduction = "integrated_umap")
ggsave(p,
       filename = paste0(results_dir, "slc17a6_gad1_featureplot.png"),
       units = "in", dpi = 600,
       height = 5, width = 10,
       bg = "white")

p <- DimPlot_scCustom(obj,
                      group.by = "DF.unadj",
                      reduction = "integrated_umap")
ggsave(p,
       filename = paste0(results_dir, "doubletfinder_umap.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

markers %>%
  filter(gene == "CLDN5") %>%
  arrange(desc(avg_log2FC))

# NEEDS INPUT: cluster ids specific to pool 149's clustering (see header).
obj@meta.data <- obj@meta.data %>%
  mutate(label1 = case_when(leiden_cluster %in% 14 ~ "Contamination",
                            leiden_cluster %in% c(6, 9, 11, 12) ~ "Inhibitory neuron",
                            leiden_cluster %in% c(1, 2, 3, 4, 5, 7, 8, 10, 13, 15, 16, 17, 18) ~ "Excitatory neuron"))

p <- DimPlot_scCustom(obj,
                      group.by = "label1",
                      reduction = "integrated_umap",
                      colors_use = c("cyan3", "yellow3", "magenta3"))
ggsave(p,
       filename = paste0(results_dir, "broad_labels.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

# Inhibitory neurons ------------------------------------------------------
message2("Subclustering inhibitory neurons")

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

p <- DimPlot_scCustom(inh,
                      reduction = "integrated_umap",
                      pt.size = 1.5)
ggsave(p,
       filename = paste0(results_dir, "inhibitory_neuron_clusters.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Transposed once up front for FindAllMarkers()'s gene-by-gene access (see
# 03_integration.R).
inh[["RNA"]]$data <- transpose_storage_order(inh[["RNA"]]$data)

inh_markers <- FindAllMarkers(inh)
write.csv(inh_markers,
          file = paste0(results_dir, "inhibitory_neuron_markers.csv"))

inh_top <- inh_markers %>%
  filter(pct.1 > 0.3 &
           str_detect(gene, "ENSG|LINC", negate = T)) %>%
  Extract_Top_Markers(5,
                      make_unique = T,
                      named_vector = F)

p <- dittoDotPlot(inh,
                  vars = inh_top,
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "inh_markers.png"),
       units = "in", dpi = 600,
       height = 5, width = 10)

{gene <- "VIP"
  p <- FeaturePlot_scCustom(inh,
                            features = gene,
                            reduction = "integrated_umap") +
    VlnPlot_scCustom(inh,
                     features = gene)
  ggsave(p,
         filename = paste0(results_dir, "inh_", gene, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 10,
         bg = "white")
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

p <- dittoBarPlot(inh,
                  var = "DF.unadj",
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "inh_doubletfinder_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 6)

# NEEDS INPUT: subcluster ids specific to pool 149's clustering.
inh@meta.data <- inh@meta.data %>%
  mutate(label2 = case_when(leiden_cluster %in% c(4, 8, 9) ~ "LAMP5 IN",
                            leiden_cluster %in% c(1, 5, 10) ~ "VIP IN",
                            leiden_cluster %in% c(2, 11) ~ "CXCL14 IN",
                            leiden_cluster %in% c(3, 12) ~ "PVALB IN",
                            leiden_cluster %in% 7 ~ "SST IN",
                            leiden_cluster == 6 ~ "Doublet"))

p <- DimPlot_scCustom(inh,
                      group.by = "label2",
                      reduction = "integrated_umap",
                      colors_use = Dark2_Pal())
ggsave(p,
       filename = paste0(results_dir, "inhibitory_neuron_labels.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Excitatory neurons ------------------------------------------------------
message2("Subclustering excitatory neurons")

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

p <- DimPlot_scCustom(ex,
                      reduction = "integrated_umap",
                      pt.size = 1.5)
ggsave(p,
       filename = paste0(results_dir, "excitatory_neuron_clusters.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

p <- dittoBarPlot(ex,
                  var = "DF.unadj",
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "ex_doubletfinder_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 8)

p <- dittoBarPlot(ex,
                  var = "tdp",
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "ex_tdp_bars.png"),
       units = "in", dpi = 600,
       height = 4, width = 8)

ex[["RNA"]]$data <- transpose_storage_order(ex[["RNA"]]$data)

ex_markers <- FindAllMarkers(ex)
write.csv(ex_markers,
          file = paste0(results_dir, "excitatory_neuron_markers.csv"))

ex_top <- ex_markers %>%
  filter(pct.1 > 0.3 &
           str_detect(gene, "ENSG|LINC", negate = T)) %>%
  Extract_Top_Markers(5,
                      make_unique = T,
                      named_vector = F)

p <- dittoDotPlot(ex,
                  vars = ex_top,
                  group.by = "leiden_cluster")
ggsave(p,
       filename = paste0(results_dir, "ex_markers.png"),
       units = "in", dpi = 600,
       height = 5, width = 10)

{gene <- "HTR2C"
  p <- FeaturePlot_scCustom(ex,
                            features = gene,
                            reduction = "integrated_umap") +
    VlnPlot_scCustom(ex,
                     features = gene)
  ggsave(p,
         filename = paste0(results_dir, "ex_", gene, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 10,
         bg = "white")
}

ex_markers %>%
  filter(cluster == 14 &
           pct.1 > 0.3 &
           str_detect(gene, "LINC|ENSG", negate = T)) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

# NEEDS INPUT: subcluster ids specific to pool 149's clustering.
ex@meta.data <- ex@meta.data %>%
  mutate(label2 = case_when(leiden_cluster %in% c(1, 2, 3, 5, 7, 8, 9, 10, 11, 14, 15) ~ "L2/3 EN",
                            leiden_cluster %in% c(6, 17) ~ "L4 EN",
                            leiden_cluster %in% c(12) ~ "L4/5 EN",
                            leiden_cluster %in% c(18, 16, 4, 13, 19) ~ "L5/6 EN"))

p <- DimPlot_scCustom(ex,
                      reduction = "integrated_umap",
                      group.by = "label2",
                      colors_use = Dark2_Pal())
ggsave(p,
       filename = paste0(results_dir, "excitatory_neuron_labels.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Fold subtype labels back and remove contamination/doublets ---------------
message2("Folding subtype labels onto the full object")

# Only label2 is carried back. Adding the subsets' whole metadata (as
# before) also overwrote the full object's leiden_cluster/seurat_clusters
# with the subclustering's ids for those cells.
new_labels <- rbind(ex@meta.data[, "label2", drop = F],
                    inh@meta.data[, "label2", drop = F])
obj <- AddMetaData(obj,
                   new_labels)

p <- DimPlot_scCustom(obj,
                      group.by = "label2",
                      reduction = "integrated_umap")
ggsave(p,
       filename = paste0(results_dir, "label2_pre-cleanup.png"),
       units = "in", dpi = 600,
       height = 5, width = 7)

# Contamination cells never get a label2 (NA), so they're dropped by the
# !is.na() check -- the old filter only dropped them implicitly.
final <- subset(obj,
                cells = colnames(obj)[!is.na(obj$label2) &
                                        obj$label2 != "Doublet"])

message2("Re-integrating the cleaned object")

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

label2_levels <- c("L2/3 EN", "L4 EN", "L4/5 EN", "L5/6 EN", "CXCL14 IN", "LAMP5 IN", "PVALB IN", "SST IN", "VIP IN")

final$label2 <- factor(final$label2,
                       levels = label2_levels)

p1 <- DimPlot_scCustom(final,
                       group.by = "label2",
                       reduction = "integrated_umap",
                       colors_use = "glasbey")
p2 <- dittoBarPlot(final,
                   var = "tdp",
                   group.by = "label2",
                   main = "Nuclear TDP content by cell type",
                   color.panel = c("navy", "orange"),
                   x.reorder = ditto_order(final$label2, label2_levels)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

p3 <- DimPlot_scCustom(final,
                       group.by = "tdp",
                       reduction = "integrated_umap",
                       pt.size = 1) +
  theme(plot.title = element_blank())

p <- p1 + p3 + p2 +
  plot_layout(ncol = 3)
ggsave(p,
       filename = paste0(results_dir, "final_labels_tdp.png"),
       units = "in", dpi = 600,
       height = 5, width = 18,
       bg = "white")

# Save ----------------------------------------------------------------------
message2("Saving barcodes, metadata, and reductions")

# Barcodes for Jackie. Cell names are now "<pool>_<sample>_<10x barcode>"
# (previously "<sample>_<10x barcode>").
write.table(colnames(final),
            file = paste0(results_dir, "filtered_barcodes.csv"),
            row.names = F,
            quote = F,
            col.names = F)

# final@meta.data %>%
#   dplyr::select(orig.ident) %>%
#   rownames_to_column(var = "barcode") %>%
#   mutate(barcode = str_split_i(barcode, "_", i = 3)) %>%
#   dplyr::rename("alias" = "orig.ident") %>%
#   write.csv()

saveRDS(final@meta.data,
        file = paste0(data_dir, "metadata.rds"))

saveRDS(final[["integrated"]],
        file = paste0(data_dir, "integrated.rds"))

saveRDS(final[["integrated_umap"]],
        file = paste0(data_dir, "integrated_umap.rds"))

# Downstream scripts rebuild the object from these on-disk pieces:
#   meta <- readRDS("data/sc_data/04_annotation/metadata.rds")
#   counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")[, rownames(meta)]
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta)
#   obj[["integrated_umap"]] <- readRDS("data/sc_data/04_annotation/integrated_umap.rds")
