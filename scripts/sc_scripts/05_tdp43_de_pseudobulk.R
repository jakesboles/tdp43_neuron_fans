# Pseudobulk DE (DESeq2 LRT and edgeR QL), TDPneg vs TDPpos, within one
# cell type, one pseudobulk sample per orig.ident. Rebuilds the cell type's
# real counts from 01_qc.R's BPCells matrix using 04_annotation.R's labels.

suppressPackageStartupMessages({
  library(DESeq2)
  library(Seurat)
  library(BPCells)
  library(tidyverse)
  library(apeglm)
  library(IHW)
  library(sva)
  library(effectsize)
  library(edgeR)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

# Cell type to test (see 05_tdp43_de_mast.R).
cell_type <- "L2/3 EN"
cell_type_file <- str_replace_all(cell_type, "[ /]", "_")

# OUTPUTS
results_dir <- paste0("results/sc_results/05_tdp43_de/", cell_type_file, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

method <- "DESeq2"

message2(paste0("Loading ", cell_type, " counts"))

meta <- readRDS("data/sc_data/04_annotation/metadata.rds")
meta <- meta[meta$label2 %in% cell_type, ]

counts <- open_matrix_dir("data/sc_data/01_qc/bpcells")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")

obj$tdp <- factor(obj$tdp,
                  levels = c("TDPpos", "TDPneg"))

# calculate cells per sample
sample_cell_nums <- obj@meta.data %>% 
  group_by(orig.ident) %>% 
  summarise(n_cells = length(orig.ident)) %>% 
  mutate(log_n_cells = log(n_cells))


# pseudobulk
bulk <- AggregateExpression(obj,
                            assays = "RNA",
                            group.by = c("orig.ident"),
                            return.seurat = F,
                            slot = "counts")

exp <- bulk$RNA # isolate counts

# format pseudobulk metadata
meta <- obj@meta.data %>%
  dplyr::select(c(orig.ident, tdp)) %>%
  distinct() %>% 
  left_join(sample_cell_nums, by = "orig.ident")
rownames(meta) <- str_replace_all(meta$orig.ident, "_", "-")

# confirm samples are aligned between counts and metadata
idx <- match(colnames(exp), rownames(meta))
meta <- meta[idx, ]

# # filter out samples with low cell counts
# log_n_cells <- meta$log_n_cells
# keep_samples <- log_n_cells > (mean(log_n_cells) - (2.5 * sd(log_n_cells)))
# exp  <- exp[, keep_samples]
# meta <- meta[keep_samples, ]

# filter for genes with at least 20 CPM in at least half of each group
sample_total_counts <- colSums(exp)
cpm_mat <- t(t(exp) / sample_total_counts) * 1e6

min_CPM <- 20
genes_keep <- rowSums(cpm_mat >= min_CPM) >= 2
# genes_keep_hc <- rowSums(cpm_mat[, hc_samples] >= min_CPM) >= (length(hc_samples) / 2)
# genes_keep_names <- intersect(names(genes_keep_pd)[genes_keep_pd == TRUE], names(genes_keep_hc)[genes_keep_hc == TRUE])
# genes_keep_lgl <- unname(genes_keep_pd & genes_keep_hc)

exp_filt <- exp[genes_keep,]

##### DIFFERENTIAL EXPRESSION #####
message2("Running differential expression")

for (method in c("DESeq2", "edgeR")){
  if (method == "DESeq2") {
    dds <- DESeqDataSetFromMatrix(countData = exp_filt,
                                  colData = meta,
                                  design = ~ tdp)
    red <- as.formula("~ 1")
    
    # fit DESeq2 model
    message("Fitting observed DESeq2 model...")
    dds_obs <- DESeq(dds,
                     test = "LRT",
                     reduced = red,
                     fitType = "local")
    
    res_obs <- results(dds_obs,
                       # contrast = c("diagnosis_2", "PD", "HC"),
                       filterFun = ihw,
                       independentFiltering = TRUE,
                       pAdjustMethod = "BH")
    
    res_df <- as.data.frame(res_obs) %>% 
      dplyr::rename(pval = pvalue,
                    padj = padj,
                    LFC = log2FoldChange)
    
    write.csv(res_df, file = paste0(results_dir, "deseq2.csv"))
    # saveRDS(dds, file = paste0(data_dir, "deseq_dataset.rds"))
    
    # Genes with FDR < 0.1, excluding genes with NA adjusted p-value
    obs_deg_count <- sum(res_df$padj < 0.05, na.rm = TRUE)
    
  } else if (method == "edgeR") {
    dgl <- DGEList(counts = exp_filt,
                   samples = meta)
    
    # library normalization
    dgl <- normLibSizes(dgl)
    
    # define model
    obs_group <- as.factor(dgl$samples$tdp)
    obs_design <- model.matrix( ~ obs_group)
    
    # estimate dispersion
    dgl_obs <- estimateDisp(dgl, 
                            design = obs_design,
                            trend.method = "locfit",
                            robust = T)
    
    # fit model and calculate effects
    fit_obs <- glmQLFit(dgl_obs, obs_design, robust = T)
    qlf_obs <- glmQLFTest(fit_obs, coef = 2)
    
    # extract and adjust p values
    pval_obs <- qlf_obs$table$PValue
    padj_obs <- p.adjust(pval_obs, method = "BH")
    qlf_obs$table$padj <- padj_obs
    
    # annotate output and save
    qlf_obs_annotated <- qlf_obs$table %>% 
      dplyr::rename(pval = PValue,
                    padj = padj,
                    LFC = logFC)
    
    write.csv(qlf_obs_annotated, file = paste0(results_dir, "edger.csv"))
    # saveRDS(dgl, file = paste0(data_dir, "dge_list_edgeR.rds"))
    
    # Genes with FDR < 0.1, excluding genes with NA adjusted p-value
    obs_deg_count <- sum(qlf_obs_annotated$padj < 0.05, na.rm = TRUE)
  }
  
  message(paste0("Genes with FDR < 0.05 from ", method, ": ", obs_deg_count))
}

# # SVA workflow
# dds <- estimateSizeFactors(dds)
# dat  <- counts(dds, normalized = TRUE)
# idx  <- rowMeans(dat) > 1
# dat  <- dat[idx, ]
# mod  <- model.matrix(~ diagnosis_2 + age + sex, colData(dds))
# mod0 <- model.matrix(~               age + sex, colData(dds))
# 
# num_SVs <- num.sv(dat, mod)
# svseq <- svaseq(dat, mod, mod0, n.sv = num_SVs)
# 
# ddssva <- dds
# included_SVs <- c()
# for (i in 1:num_SVs) {
#   SV_batch_anova <- aov(as.formula(paste0("svseq$sv[,", i, "] ~ colData(ddssva)$batch_id")))
#   SV_batch_assoc_p <- summary(SV_batch_anova)[[1]]$`Pr(>F)`[1]
#   
#   if (SV_batch_assoc_p < 0.1) {
#     ddssva[[paste0("SV", i)]] <- svseq$sv[,i]
#     included_SVs <- append(included_SVs, i)
#   }
# }
# 
# pdf(paste0(output_dir, "qc/", celltype, "_SVA_plots.pdf"))
# for (i in 1:num_SVs) {
#   stripchart(svseq$sv[, i] ~ dds$batch_id, vertical = TRUE, main = paste0("SV", i, ifelse(i %in% included_SVs, " *", "")))
#   abline(h = 0)
# }
# dev.off()
# 
# design(ddssva) <- as.formula(paste0("~ diagnosis_2 + age + sex + ", 
#                                     paste(paste0("SV", included_SVs), collapse = " + ")))
# red <- as.formula(paste0("~ 1 + age + sex + ", 
#                          paste(paste0("SV", included_SVs), collapse = " + ")))
# 
# ddssva <- DESeq(ddssva,
#                 test = "LRT",
#                 reduced = red,
#                 fitType = "local")
# 
# # residual <- mcols(dds)$dispGeneEst - mcols(dds)$dispFit
# # mad(residual)
# # median(abs(residual))
# 
# res <- results(ddssva,
#                contrast = c("diagnosis_2", "PD", "HC"),
#                filterFun = ihw,
#                independentFiltering = T,
#                pAdjustMethod = "BH")
# 
# # sum(res$padj < 0.1, na.rm=TRUE)
# # sum(res$padj < 0.05, na.rm=TRUE)
# 
# res_df <- as.data.frame(res)
# 
# write.csv(res_df, 
#           file = paste0(output_dir, celltype, "_deseq.csv"))
# 
# # res_lfc_shrink <- lfcShrink(dds, coef = "diagnosis_2_PD_vs_HC",
# #                             type = "apeglm")
# # 
# # res_lfc_shrink_df <- as.data.frame(res_lfc_shrink)
# # write.csv(res_lfc_shrink_df,
# #           file = paste0(output_dir, celltype, "_deseq_lfc_shrunk.csv"))
# 
# # PCA plot
# vsd <- vst(dds, nsub = min(1000, dim(dds)[1]), blind = TRUE)
# 
# pca_df <- plotPCA(
#   vsd,
#   intgroup = c("diagnosis_2", "batch_id", "age", "sex"),
#   returnData = TRUE
# ) %>% 
#   rownames_to_column("sample_id") %>% 
#   left_join(sample_cell_nums, by = "sample_id")
# 
# write_csv(pca_df, paste0(output_dir, "qc/", celltype, "_pca_df.csv"))
# 
# percentVar <- round(100 * attr(pca_df, "percentVar"))
# 
# pdf(paste0(output_dir, "qc/", celltype, "_PCA_plot.pdf"))
# ggplot(pca_df, aes(PC1, PC2, color = batch_id, shape = diagnosis_2)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   labs(
#     x = paste0("PC1: ", percentVar[1], "%"),
#     y = paste0("PC2: ", percentVar[2], "%")
#   )
# dev.off()
# 
# # correlation grid for PCs
# pc_meta_stats <- function(pc, meta) {
#   if (is.numeric(meta)) {
#     ct <- cor.test(pc, meta)
#     c(assoc = unname(ct$estimate), p = ct$p.value)
#     
#   } else if (is.factor(meta) && nlevels(meta) == 2) {
#     meta_num <- as.numeric(meta) - 1
#     ct <- cor.test(pc, meta_num)
#     c(assoc = unname(ct$estimate), p = ct$p.value)
#     
#   } else if (is.factor(meta)) {
#     a <- aov(pc ~ meta)
#     eta2 <- eta_squared(a)
#     c(assoc = sqrt(eta2$Eta2), p = summary(a)[[1]]$`Pr(>F)`[1])
#     
#   } else {
#     c(assoc = NA_real_, p = NA_real_)
#   }
# }
# 
# pcs <- c("PC1", "PC2")
# meta_cols <- c("batch_id", "age", "sex", "n_cells", "diagnosis_2")
# 
# stats <- expand.grid(
#   metadata = meta_cols,
#   PC = pcs,
#   stringsAsFactors = FALSE
# )
# 
# stats <- stats |>
#   dplyr::rowwise() |>
#   dplyr::mutate(
#     res = list(
#       pc_meta_stats(
#         pca_df[[PC]],
#         pca_df[[metadata]]
#       )
#     ),
#     association = res["assoc"],
#     p_value = res["p"]
#   ) |>
#   dplyr::ungroup() |>
#   dplyr::select(-res)
# 
# stats$padj <- p.adjust(stats$p_value, method = "BH")
# 
# stats$sig <- cut(
#   stats$padj,
#   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
#   labels = c("***", "**", "*", ""),
#   right = FALSE
# )
# 
# stats$label <- ifelse(
#   is.na(stats$association),
#   "",
#   paste0(
#     sprintf("%.3f", stats$association),
#     stats$sig
#   )
# )
# 
# pdf(paste0(output_dir, "qc/", celltype, "_PC_correlations.pdf"))
# ggplot(stats, aes(PC, metadata, fill = association)) +
#   geom_tile(color = "white") +
#   geom_text(
#     aes(label = label),
#     size = 4
#   ) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     na.value = "grey90",
#     name = "Effect size"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(size = 10),
#     panel.grid = element_blank()
#   )
# dev.off()



# ### Test
# gene_res <- tibble(gene_id = rownames(res_df),
#                    padj = res_df$padj) %>% 
#   mutate(sig = padj < 0.1) %>% 
#   select(-padj)
# 
# gene_res_old <- tibble(gene_id = res_df_old$...1,
#                        padj = res_df_old$padj) %>% 
#   mutate(sig = padj < 0.1) %>% 
#   select(-padj)
# 
# compare <- rbind(gene_res, gene_res_old) %>% 
#   pivot_wider(id_cols = )
# 
# # R example
# merged <- merge(gene_res, gene_res_old, by = "gene_id", suffixes = c("_A", "_B"))
# 
# table_sig <- table(merged$sig_A, merged$sig_B)
# fisher.test(table_sig)


# hist(res_df$pvalue, breaks = 20)
# hist(res_df$padj, breaks = 20)
# 
# plotCounts(dds, gene = "SLC39A8",
#            intgroup = "diagnosis_2")
# 
# obj <- NormalizeData(obj)
# 
# VlnPlot(obj,
#         group.by = "diagnosis_2",
#         features = "SLC39A8")
# 
# Idents(obj) <- "diagnosis_2"
# 
# de <- FindMarkers(obj,
#                   ident.1 = "PD",
#                   ident.2 = "HC")
# 
# de["SLC39A8",]
# 
# de %>%
#   filter(pct.1 > 0.3) %>%
#   arrange(desc(avg_log2FC)) %>%
#   head(30)
# 
# de %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val_adj))) + 
#   geom_point()
# 
# res_df %>% 
#   ggplot(aes(x = log2FoldChange,
#              y = -log10(padj))) + 
#   geom_point() + 
#   geom_hline(yintercept = -log10(0.05))
