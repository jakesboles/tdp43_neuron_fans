# DESeq2 on the bulk RNA-seq of TDP-43- vs TDP-43+ neuronal nuclei from
# GEO GSE126542 (Liu et al.), for comparison with this dataset's DE.

suppressMessages({
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(IHW)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/tdp43_neuron_fans")

results_dir <- "results/sc_results/06_liu_de/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

# NEEDS INPUT: external_data/ is gitignored -- make sure this GEO counts
# file is in external_data/ under the project root.
message2("Reading GSE126542 counts")

df <- read_table("external_data/GSE126542_NeuronalNuclei_RNAseq_counts.txt")

genes <- df$S1_TN

col.names <- read_table("external_data/GSE126542_NeuronalNuclei_RNAseq_counts.txt",
                        n_max = 1)
col.names <- colnames(col.names)

df <- read_table("external_data/GSE126542_NeuronalNuclei_RNAseq_counts.txt",
                 skip = 1)

df <- rbind(colnames(df), df)
colnames(df) <- c("gene", col.names)

df <- df %>%
  na.omit()

df <- df %>%
  column_to_rownames(var = "gene")

counts <- data.matrix(df)

meta <- data.frame(
  sample = col.names
) %>%
  mutate(donor = str_split_i(sample, "_", i = 1),
         tdp = str_split_i(sample, "_", i = 2) %>%
           factor(levels = c("TN", "TP"),
                  labels = c("TDPneg", "TDPpos")))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design= ~ tdp)

dds

keep <- rowSums(counts(dds) >= 10) >= 7
dds <- dds[keep,]

dds$tdp

message2("Running DESeq2")

dds <- DESeq(dds)

res <- results(dds,
               independentFiltering = T,
               filterFun = ihw)
head(as.data.frame(res))

res %>%
  as.data.frame() %>% 
  write.csv(paste0(results_dir, "deseq.csv"))
