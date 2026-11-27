# Builds one Seurat object from the Cell Ranger per-sample outputs of every
# 10x pool listed in config/pools.tsv, computes QC metrics, plots them, and
# applies the QC filter. Counts are saved on-disk with BPCells, plus the
# metadata as RDS, for every later script to rebuild the object from.

suppressMessages({
  library(Seurat)
  library(BPCells)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(ggbeeswarm)
  library(janitor)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

results_dir <- "results/sc_results/01_qc/"
data_dir <- "data/sc_data/01_qc/"

for (i in c(results_dir, data_dir)){
  dir.create(i,
             showWarnings = F,
             recursive = T)
}

# QC thresholds -------------------------------------------------------------
max_mito <- 10
min_count <- 1500
max_count <- 80000
min_feature <- 1000
max_feature <- 8500

# Find Cell Ranger outputs for every pool -----------------------------------
message2("Finding Cell Ranger outputs")

# Pools without short-read data yet have cellranger_dir = NA and are skipped.
pools <- read.delim("config/pools.tsv",
                    na.strings = "NA") %>%
  filter(!is.na(cellranger_dir))

if (nrow(pools) == 0){
  stop("No pool in config/pools.tsv has a cellranger_dir set.")
}

samples <- map(seq_len(nrow(pools)), function(i){
  dirs <- list.dirs(pools$cellranger_dir[i],
                    recursive = F,
                    full.names = T)
  if (length(dirs) == 0){
    stop(paste0("No per-sample folders found in ", pools$cellranger_dir[i],
                " (pool ", pools$pool[i], ")"))
  }
  data.frame(pool = pools$pool[i],
             sample = basename(dirs),
             path = paste0(dirs, "/sample_filtered_feature_bc_matrix/"))
}) %>%
  list_rbind() %>%
  # Pool-prefixed so sample ids stay unique when two pools reuse a
  # sample name (e.g. the same donor sorted in two pools).
  mutate(id = paste0(pool, "_", sample))

missing <- samples$path[!dir.exists(samples$path)]
if (length(missing) > 0){
  stop(paste0("Missing Cell Ranger matrix folder(s):\n",
              paste(missing, collapse = "\n")))
}

print(samples)

# Create Seurat objects -----------------------------------------------------
message2("Creating Seurat objects")

objs <- list()

for (i in seq_len(nrow(samples))){

  message(samples$id[i])

  counts <- Read10X(samples$path[i])
  objs[[i]] <- CreateSeuratObject(counts = counts,
                                  project = samples$id[i])
  objs[[i]]$pool <- samples$pool[i]
  objs[[i]]$sample <- samples$sample[i]
}

message2("Merging Seurat objects")

# Cell names become "<pool>_<sample>_<10x barcode>".
obj <- Merge_Seurat_List(objs,
                         add.cell.ids = samples$id)

obj@meta.data %>%
  group_by(pool, orig.ident) %>%
  summarize(n = n(),
            .groups = "drop") %>%
  janitor::adorn_totals() %>%
  write.csv(paste0(results_dir, "cell_counts.csv"),
            quote = F,
            row.names = F)

obj <- JoinLayers(obj)

# Cell Ranger sample names are "<tissue>_<tdp>", where tissue is the brain
# region (e.g. Frontal_TDPneg; see config/jsb149_cellranger_config.csv).
obj@meta.data <- obj@meta.data %>%
  mutate(tissue = str_split_i(sample, "_", i = 1),
         tdp = str_split_i(sample, "_", i = 2))

obj <- obj %>%
  Add_Mito_Ribo("human") %>%
  Add_Cell_Complexity()

# Raw QC plots --------------------------------------------------------------
message2("Plotting raw QC metrics")

raw_plot <- function(feature, log = F){
  p <- VlnPlot_scCustom(obj,
                        features = feature,
                        group.by = "tissue",
                        split.by = "tdp") +
    theme(axis.title.x = element_blank())

  if (isTRUE(log)){
    p <- p +
      scale_y_log10()
  }

  return(p)
}

p <- raw_plot("nCount_RNA", log = T)
ggsave(p,
       filename = paste0(results_dir, "raw_ncount.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

p <- raw_plot("nFeature_RNA", log = T)
ggsave(p,
       filename = paste0(results_dir, "raw_nfeature.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

p <- raw_plot("percent_mito")
ggsave(p,
       filename = paste0(results_dir, "raw_mito.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

p <- raw_plot("log10GenesPerUMI")
ggsave(p,
       filename = paste0(results_dir, "raw_complexity.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Lines drawn at the nFeature thresholds actually applied below (this plot
# previously drew 9000 while the filter used 8500).
p <- VlnPlot(obj,
             features = "nFeature_RNA",
             log = T) +
  geom_hline(yintercept = c(min_feature, max_feature))
ggsave(p,
       filename = paste0(results_dir, "raw_nfeature_thresholds.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Apply QC filter -----------------------------------------------------------
message2("Applying QC filter")

obj@meta.data <- obj@meta.data %>%
  mutate(discard = if_else(percent_mito > max_mito |
                             nCount_RNA < min_count |
                             nCount_RNA > max_count |
                             nFeature_RNA < min_feature |
                             nFeature_RNA > max_feature,
                           T, F))

table(obj$orig.ident, obj$discard)

obj@meta.data %>%
  janitor::tabyl(orig.ident, discard) %>%
  write.csv(paste0(results_dir, "discard_counts.csv"),
            quote = F,
            row.names = F)

meta <- obj@meta.data %>%
  mutate(discard = factor(discard,
                          levels = c(F, T),
                          labels = c("Keep", "Discard")))

# NEEDS INPUT: with levels Keep, Discard this colors kept cells red and
# discarded cells black -- swap the two colors if that's backwards.
discard_plot <- function(feature,
                         log = F){
  p <- meta %>%
    ggplot(aes(x = orig.ident,
               y = {{feature}})) +
    geom_quasirandom(aes(color = discard)) +
    scale_color_manual(values = c("firebrick1", "black")) +
    ggtitle(deparse(substitute(feature))) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if (isTRUE(log)){
    p <- p +
      scale_y_log10()
  }

  return(p)

}

# Widen with the number of samples so the x axis stays readable.
discard_width <- max(5, 1.25 * length(unique(meta$orig.ident)))

p <- discard_plot(nCount_RNA,
             log = T)
ggsave(p,
       filename = paste0(results_dir, "filtered_ncount.png"),
       units = "in", dpi = 600,
       height = 4, width = discard_width)

p <- discard_plot(nFeature_RNA,
             log = T)
ggsave(p,
       filename = paste0(results_dir, "filtered_nfeature.png"),
       units = "in", dpi = 600,
       height = 4, width = discard_width)

p <- discard_plot(percent_mito,
             log = F)
ggsave(p,
       filename = paste0(results_dir, "filtered_mito.png"),
       units = "in", dpi = 600,
       height = 4, width = discard_width)

p <- discard_plot(log10GenesPerUMI)
ggsave(p,
       filename = paste0(results_dir, "filtered_complexity.png"),
       units = "in", dpi = 600,
       height = 4, width = discard_width)

obj <- subset(obj,
              subset = discard == F)

# Save ----------------------------------------------------------------------
message2("Saving filtered counts as BPCells on-disk matrix")

# Real raw counts for the QC-filtered cells. Every later script rebuilds its
# object from this matrix (subset to its own cells) plus its own metadata,
# so normalized data never has to be written to disk.
bpcells_dir <- paste0(data_dir, "bpcells")
if (dir.exists(bpcells_dir)){
  unlink(bpcells_dir, recursive = T)
}

# dgCMatrix stores counts as doubles; write_matrix_dir() needs an explicitly
# integer-typed matrix.
counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")

write_matrix_dir(mat = counts_out,
                 dir = bpcells_dir)

message2("Saving metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_dir, "metadata.rds"))

# Downstream scripts rebuild the object from these on-disk pieces:
#   counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
#   meta <- readRDS("data/sc_data/01_qc/metadata.rds")
#   obj <- CreateSeuratObject(counts = counts[, rownames(meta)], meta.data = meta)
