library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(janitor)

setwd("/projects/b1169/boles/JSB149")

plots_dir <- "plots/01_qc/"
tab_dir <- "tab_data/01_qc/"
data_dir <- "data/01_qc/"

for (i in c(plots_dir, tab_dir, data_dir)){
  dir.create(i,
             showWarnings = F,
             recursive= T)
}

cellranger <- list.dirs("/projects/b1042/Gate_Lab/boles/JSB149/cellranger/JSB149/outs/per_sample_outs",
                        recursive = F,
                        full.names = T)

samples <- str_split_i(cellranger, "/", i = 11)

cellranger <- paste0(cellranger, "/sample_filtered_feature_bc_matrix/")

objs <- list()

for (i in seq_along(cellranger)){
  
  message(samples[i])
  
  counts <- Read10X(cellranger[i])
  objs[[i]] <- CreateSeuratObject(counts = counts,
                                  project = samples[i])
}

obj <- Merge_Seurat_List(objs,
                         add.cell.ids = samples)

obj@meta.data %>% 
  group_by(orig.ident) %>% 
  summarize(n = n()) %>%
  janitor::adorn_totals() %>%
  write.csv(paste0(tab_dir, "cell_counts.csv"),
            quote = F,
            row.names = F)
  

obj <- JoinLayers(obj)

obj@meta.data <- obj@meta.data %>% 
  mutate(tissue = str_split_i(orig.ident, "_", i = 1),
         tdp = str_split_i(orig.ident, "_", i = 2))

obj <- obj %>%
  Add_Mito_Ribo("human") %>%
  Add_Cell_Complexity()

VlnPlot_scCustom(obj,
                 features = "nCount_RNA",
                 group.by = "tissue",
                 split.by = "tdp") + 
  scale_y_log10() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plots_dir, "raw_ncount.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

VlnPlot_scCustom(obj,
                 features = "nFeature_RNA",
                 group.by = "tissue",
                 split.by = "tdp") + 
  scale_y_log10() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plots_dir, "raw_nfeature.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

VlnPlot_scCustom(obj,
                 features = "percent_mito",
                 group.by = "tissue",
                 split.by = "tdp") + 
  # scale_y_log10() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plots_dir, "raw_mito.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

VlnPlot_scCustom(obj,
                 features = "log10GenesPerUMI",
                 group.by = "tissue",
                 split.by = "tdp") + 
  # scale_y_log10() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plots_dir, "raw_complexity.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

VlnPlot(obj,
        features = "nFeature_RNA",
        log = T) + 
  geom_hline(yintercept = 9000)


obj@meta.data <- obj@meta.data %>% 
  mutate(discard = if_else(percent_mito > 10 | 
                             nCount_RNA < 1500 | 
                             nCount_RNA > 80000 | 
                             nFeature_RNA < 1000 | 
                             nFeature_RNA > 8500,
                           T, F))

table(obj$orig.ident, obj$discard)

meta <- obj@meta.data %>%
  mutate(discard = factor(discard,
                          levels = c(F, T),
                          labels = c("Keep", "Discard")))

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

discard_plot(nCount_RNA,
             log = T)
ggsave(filename = paste0(plots_dir, "filtered_ncount.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

discard_plot(nFeature_RNA,
             log = T)
ggsave(filename = paste0(plots_dir, "filtered_nfeature.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

discard_plot(percent_mito,
             log = F)
ggsave(filename = paste0(plots_dir, "filtered_mito.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

discard_plot(log10GenesPerUMI)
ggsave(filename = paste0(plots_dir, "filtered_complexity.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

obj <- subset(obj,
              subset = discard == F)

obj@meta.data %>%
  group_by(orig.ident) %>%
  summarize(n = n()) %>%
  write.table("/projects/p31535/boles/JSB149/02_doubletfinder_params.txt",
              quote = F,
              sep = ",",
              row.names = F,
              col.names = F)

saveRDS(obj,
        file = paste0(data_dir, "obj.rds"))
