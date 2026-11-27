# Runs DoubletFinder (adjusted and unadjusted) on one sample, loading that
# sample's QC-filtered counts from 01_qc.R's BPCells matrix. Runs as a SLURM
# job array (see jobs/sc_jobs/02_doubletfinder.sh), one task per sample.
# Script again adapted from Anne's projects.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(Matrix)
  library(DoubletFinder)
  library(scCustomize)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

results_dir <- "results/sc_results/02_doubletfinder/"
dir.create(results_dir,
           showWarnings = F, recursive = T)

# Figure out which sample this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/sc_jobs/02_doubletfinder.sh), one task per ",
       "sample, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

message2("Loading metadata")
meta_all <- readRDS("data/sc_data/01_qc/metadata.rds")

# Sorting guarantees the same sample <-> array index mapping every run.
samples <- sort(unique(meta_all$orig.ident))

if (task_id < 1 | task_id > length(samples)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(samples), " samples -- check the --array range in ",
              "jobs/sc_jobs/02_doubletfinder.sh against ",
              "length(unique(orig.ident))."))
}

sample_id <- samples[task_id]

message2(paste0("Processing ", sample_id, " (task ", task_id, "/",
                length(samples), ")"))

# Load this sample's counts -------------------------------------------------

meta_sample <- meta_all[meta_all$orig.ident == sample_id, ]

mat <- open_matrix_dir("data/sc_data/01_qc/bpcells")
mat <- mat[, rownames(meta_sample)]

# DoubletFinder's paramSweep() samples cell pairs WITH replacement, and
# BPCells' `[` rejects duplicate column selections, so materialize this one
# sample in memory first.
mat <- as(mat, "dgCMatrix")

obj <- CreateSeuratObject(counts = mat, meta.data = meta_sample)

# NEEDS INPUT: this is the post-QC cell count, as in the original script
# (it came from 01_qc.R's params file written after filtering). Doublets
# form at loading, so the pre-QC count -- or, with OCM, possibly the whole
# channel's cell count -- may be the better basis. 0.04 per 10,000 cells
# matches 10x's GEM-X rate.
cell_count <- ncol(obj)
doublet_rate <- (cell_count / 10000) * 0.04

# Define function to run DoubletFinder ---------------------------------------

run_doubletfinder <- function(s, doublet_rate) {

  # Standard normalization and scaling
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Default Seurat clustering and UMAP
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunUMAP(s, dims = 1:10) # had to set check_duplicates to FALSE for some samples

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # Run DoubletFinder without homotypic adjustment
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"

  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  print(pANN)
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)

  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"

  p1 <- DimPlot_scCustom(s,
                         reduction = "umap",
                         group.by = "DF.unadj",
                         pt.size = 1) +
    ggtitle("Unadjusted")

  p2 <- DimPlot_scCustom(s,
                         reduction = "umap",
                         group.by = "DF.adj",
                         pt.size = 1) +
    ggtitle("Adjusted for homotypic proportion")

  p <- p1 + p2

  ggsave(p,
         filename = paste0(results_dir, sample_id, "_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 10)

  return(s)
}

message2("Running DoubletFinder")

obj <- run_doubletfinder(obj, doublet_rate)

# Save ----------------------------------------------------------------------
message2("Saving DoubletFinder calls")

# 03_integration.R reads these per-sample tables and checks every sample is
# present. The per-sample Seurat object is no longer saved: nothing
# downstream read it, and the counts already live in 01_qc.R's matrix.
write.csv(obj@meta.data,
          file = paste0(results_dir, sample_id, "_meta.csv"))
