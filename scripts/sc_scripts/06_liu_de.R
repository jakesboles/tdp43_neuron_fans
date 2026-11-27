library(tidyverse)
library(DESeq2)
library(apeglm)
library(IHW)

setwd("/projects/b1169/boles/JSB149")

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

dds <- DESeq(dds)

res <- results(dds,
               independentFiltering = T,
               filterFun = ihw)
View(as.data.frame(res))

dir.create("tab_data/06_liu_de",
           showWarnings = F,
           recursive = T)

res %>%
  as.data.frame() %>% 
  write.csv("tab_data/06_liu_de/deseq.csv")
