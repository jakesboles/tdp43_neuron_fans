# MAST single-cell DE (TDPneg vs TDPpos) within one cell type, with cellular
# detection rate as a covariate and a per-sample random intercept. Rebuilds
# the cell type's counts from 01_qc.R's BPCells matrix using
# 04_annotation.R's labels.

suppressMessages({
  library(Seurat)
  library(BPCells)
  library(MAST)
  library(tidyverse)
  library(patchwork)
  library(lme4)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

p_thresh <- 0.05
fc_thresh <- 0.25

# Cell type to test; also set in the other 05_tdp43_de_* scripts and the
# scripts that read their results (07, summary_plots, Novartis comparison).
cell_type <- "L2/3 EN"
cell_type_file <- str_replace_all(cell_type, "[ /]", "_")

# Shared with 05_tdp43_de_pseudobulk.R and 05_tdp43_de_wilcox.R
output_dir <- paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load expression within cell type
message2(paste0("Loading ", cell_type, " counts"))

meta <- readRDS("data/sc_data/04_annotation/metadata.rds")
meta <- meta[meta$label2 %in% cell_type, ]

counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

# One cell type fits in memory, and MAST needs a plain matrix anyway.
counts <- as(counts, "dgCMatrix")

s <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
s <- NormalizeData(s)


#### Differential Expression: Single Cell ####
groupvar = "tdp"
cond_interest = "TDPneg"
cond_reference = "TDPpos"

Idents(s) <- groupvar

if (!(cond_interest %in% unique(Idents(s)) & cond_reference %in% unique(Idents(s)))) {
  stop("Both conditions not present")
}

# Calculate cellular detection rate (CDR) as MAST covariate
cdr <- colMeans(GetAssayData(s, assay = "RNA", layer = "data") > 0)
s$cdr <- cdr
s$cdr_centered <- scale(s$cdr)

# Get avg_log2FC and percent expression for all genes
LFC <- FoldChange(object = s, 
                  ident.1 = cond_interest, ident.2 = cond_reference, 
                  fc.name = "avg_log2FC", base = 2,
                  assay = "RNA", layer = "data") %>% data.frame()

# Calculate overall percent expression
overall <- rowMeans(GetAssayData(s, assay = "RNA", layer = "data") > 0)
overall <- data.frame(gene = names(overall), pct = unlist(overall))
print(sum(rownames(overall) != rownames(LFC)))

# Gene expression filtering
# Either it's in >1% overall and >5% of one group OR it's >5% of both groups
genes_keep <- row.names(LFC)[((LFC$pct.1 >= 0.05 | LFC$pct.2 >= 0.05) & overall$pct >= 0.01) | (LFC$pct.1 >= 0.05 & LFC$pct.2 >= 0.05)]

# Remove contamination
if (length(grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)) != 0) {
  genes_keep <- genes_keep[-grep(pattern = "^RPS|^RPL|^MT-|^HB", x = genes_keep)] 
}

print(paste0("Testing ", length(genes_keep), " genes"))

# Get all expression data and filter
expressionmat_full <- GetAssayData(s, assay = "RNA", layer = 'data')
expressionmat <- expressionmat_full[genes_keep,]
expressionmat <- as.matrix(expressionmat)

# Cell-level and feature-level meta data for MAST
cdat <- s@meta.data
cdat$wellKey <- row.names(cdat)
fdat <- data.frame(primerid = row.names(expressionmat), row.names = row.names(expressionmat))

# Make sure columns of expression data match rows of meta data
print(sum(colnames(expressionmat) != row.names(cdat)))

# Create SingleCellAssay object
sca <- FromMatrix(expressionmat, cdat, fdat, check_sanity = FALSE)

# Set reference level for comparison 
group_ordered <- factor(colData(sca)[,groupvar])
group_ordered <- relevel(group_ordered, cond_reference)
colData(sca)$group <- group_ordered

# Enable parallel processing 
options(mc.cores = 16)

# MAST does not allow ebayes (empirical Bayes estimation to regularize variance) for glmer
message2("Fitting MAST model")

zlm_group <- zlm(~ tdp + cdr_centered + (1 | orig.ident), 
                 sca, 
                 ebayes = FALSE, 
                 method = "glmer", 
                 # hook = continuous_residuals_hook,
                 parallel = T, 
                 fitArgsD = list(nAGQ = 0))

lrt_name <- paste0(groupvar, cond_reference)
summary <- summary(zlm_group, doLRT = lrt_name, parallel = T) # FindMarkers: doLRT = "conditionGroup2"

summary_data <- summary$datatable %>% data.frame()
p_val <- summary_data[summary_data[, "component"] == "H", 4]
genes.return <- summary_data[summary_data[, "component"] == "H", 1]

# Compile results
results <- data.frame("p_val" = p_val, "gene" = genes.return, 
                      row.names = genes.return)
results$BH <- p.adjust(results$p_val, method = "BH")

LFC_use <- LFC[match(results$gene, row.names(LFC)),]
print(sum(row.names(LFC_use) != results$gene))
results$avg_log2FC <- LFC_use$avg_log2FC

results$DE <- "Not DE"
results$DE[results$avg_log2FC > fc_thresh & results$BH < p_thresh] <- "Upregulated"
results$DE[results$avg_log2FC < -fc_thresh & results$BH < p_thresh] <- "Downregulated"

# For volcano plot labels
results$DE_gene <- NA
results$DE_gene[results$DE != "Not DE"] <- results$gene[results$DE != "Not DE"]

table(results$DE)

# Save results
message2("Saving results")

write.csv(results, paste0(output_dir, "mast.csv"))
