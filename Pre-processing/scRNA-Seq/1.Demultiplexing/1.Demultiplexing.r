---
#title: "Redmond et al. 2021 - Step 1: Demultiplexing"
---

# Setup

library(Seurat)
library(ShortRead)
library(deMULTIplex)
library(KernSmooth)
library(reshape2)
library(Rtsne)
library(stringdist)
library(knitr)
library(dplyr)
library(ggplot2)
library(deMULTIplex)
library(future)


# Reading the Cellranger output and Creating a Seurat Object

svz_matrix <- Read10X(
  data.dir = c("/data4/svz_10x/wholecell/data/pipestance/g1/outs/filtered_feature_bc_matrix",
               "/data4/svz_10x/wholecell/data/pipestance/g2/outs/filtered_feature_bc_matrix"),
  unique.features = TRUE,
  strip.suffix = TRUE
)

exp <- CreateSeuratObject(counts = svz_matrix$`Gene Expression`)
exp[['Protein']] <- CreateAssayObject(counts = svz_matrix$`Antibody Capture`)
saveRDS(exp, file = "experiment_raw.rds")


# 1. Demultiplexing

#Determining the cell ids in 10x chip lane 1 and lane 2
# Multi-seq barcodes used: A3(TGAGACCT), A4(GCACACGC), A5(AGAGAGAG), A6 (TCACAGCA)
multi.ref <- c("TGAGACCT", "GCACACGC", "AGAGAGAG", "TCACAGCA")

all.cell.ids <- names(exp$orig.ident)

#A typical cell id in lane 1 looks like this: "ATTCCTACATACTGAC" (16 bases). However, cell ids in lane 2 have a "2_" prefix, looking like this: "2_AAACCCAAGAGCTGAC". This prefix is added on cellranger and we have to remove it in order to perform the alignment of barcode reads.

b1.cell.ids <- all.cell.ids[!grepl("^2_", all.cell.ids)]
b2.cell.ids <- all.cell.ids[grepl("^2_", all.cell.ids)]

exp@meta.data[b1.cell.ids,"Lane"] <- "Lane_1"
exp@meta.data[b2.cell.ids,"Lane"] <- "Lane_2"

#we have to remove the 2_ prefix of b2 ids in order to align reads:
b2.cell.ids <- sub("2_", "", b2.cell.ids)


## Loading FastQ files
#tgacaat - Barcode library lane 1
#acagtgat - Barcode library lane2

# B1:
readTable.b1 <- MULTIseq.preProcess(R1 = "../fastq.files/TGACCAAT_S5_L001_R1_001_subsampled.fastq.gz", R2 = "../fastq.files/TGACCAAT_S5_L001_R2_001_subsampled.fastq.gz", cellIDs = b1.cell.ids)
saveRDS(readTable.b1, file= "readTable.b1.rds")

# B2:
readTable.b2 <- MULTIseq.preProcess(R1 = "../fastq.files/ACAGTGAT_S6_L001_R1_001_subsampled.fastq.gz", R2 = "../fastq.files/ACAGTGAT_S6_L001_R2_001_subsampled.fastq.gz", cellIDs = b2.cell.ids)
saveRDS(readTable.b2, file= "readTable.b2.rds")


## Perform MULTI-seq sample barcode alignment
bar.table.b1 <- MULTIseq.align(readTable.b1, b1.cell.ids, multi.ref)
saveRDS(bar.table.b1, file = "bar.table.b1.rds")

bar.table.b2 <- MULTIseq.align(readTable.b2, b2.cell.ids, multi.ref)
saveRDS(bar.table.b2, file = "bar.table.b2.rds")


## Visual inspection in the barcode space
#B1
bar.tsne.b1 <- barTSNE(bar.table.b1 [,1:4]) 
saveRDS(bar.tsne.b1, file="bar.tsne.b1.rds")


for (i in 3:ncol(bar.tsne.b1)) {
    g <- ggplot(bar.tsne.b1, aes(x = TSNE1, y = TSNE2, color = bar.tsne.b1[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne.b1)[i]) +
    theme(legend.position = "none")  + coord_fixed()
    print(g)
}
dev.off()

#B2
bar.tsne.b2 <- barTSNE(bar.table.b2 [,1:4]) 
saveRDS(bar.tsne.b2, file="bar.tsne.b2.rds")

for (i in 3:ncol(bar.tsne.b2)) {
    g <- ggplot(bar.tsne.b2, aes(x = TSNE1, y = TSNE2, color = bar.tsne.b2[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne.b2)[i]) +
    theme(legend.position = "none")  + coord_fixed()
    print(g)
}
dev.off()


## Round 1 / B1
## Perform Quantile Sweep
bar.table.b1.full <- bar.table.b1[,1:4]
good.bars <- paste("Bar",1:4,sep="")  # Selecting all 4 barcodes
bar.table.b1 <- bar.table.b1.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list.b1 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b1[[n]] <- classifyCells(bar.table.b1, q=q)
  names(bar.table_sweep.list.b1)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results.b1.1 <- findThresh(call.list=bar.table_sweep.list.b1)
ggplot(data=threshold.results.b1.1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b1.1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

# Finalize round 1 classifications, remove negative cells
b1.round1.calls <- classifyCells(bar.table.b1, q=findQ(threshold.results.b1.1$res, threshold.results.b1.1$extrema))
b1.neg.cells <- names(b1.round1.calls)[which(b1.round1.calls == "Negative")]
bar.table.b1 <- bar.table.b1[-which(rownames(bar.table.b1) %in% b1.neg.cells), ]


## Round 1 / B2 
## Perform Quantile Sweep
bar.table.b2.full <- bar.table.b2[,1:4]
good.bars <- paste("Bar",1:4,sep="")  # Selecting all 4 barcodes
bar.table.b2 <- bar.table.b2.full[, good.bars]  # Remove missing bars and summary columns
bar.table_sweep.list.b2 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b2[[n]] <- classifyCells(bar.table.b2, q=q)
  names(bar.table_sweep.list.b2)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results.b2.1 <- findThresh(call.list=bar.table_sweep.list.b2)
ggplot(data=threshold.results.b2.1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b2.1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

# Finalize round 1 classifications, remove negative cells
b2.round1.calls <- classifyCells(bar.table.b2, q=findQ(threshold.results.b2.1$res, threshold.results.b2.1$extrema))
b2.neg.cells <- names(b2.round1.calls)[which(b2.round1.calls == "Negative")]
bar.table.b2 <- bar.table.b2[-which(rownames(bar.table.b2) %in% b2.neg.cells), ]


## Round 2 /B1

bar.table_sweep.list.b1 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b1[[n]] <- classifyCells(bar.table.b1, q=q)
  names(bar.table_sweep.list.b1)[n] <- paste("q=",q,sep="")
}

threshold.results.b1.2 <- findThresh(call.list=bar.table_sweep.list.b1)
b1.round2.calls <- classifyCells(bar.table.b1, q=findQ(threshold.results.b1.2$res, threshold.results.b1.2$extrema))
b1.neg.cells <- c(b1.neg.cells, names(b1.round2.calls)[which(b1.round2.calls == "Negative")])

ggplot(data=threshold.results.b1.2$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b1.2$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))


## Round 2 /B2
bar.table_sweep.list.b2 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b2[[n]] <- classifyCells(bar.table.b2, q=q)
  names(bar.table_sweep.list.b2)[n] <- paste("q=",q,sep="")
}

threshold.results.b2.2 <- findThresh(call.list=bar.table_sweep.list.b2)
b2.round2.calls <- classifyCells(bar.table.b2, q=findQ(threshold.results.b2.2$res, threshold.results.b2.2$extrema))
b2.neg.cells <- c(b2.neg.cells, names(b2.round2.calls)[which(b2.round2.calls == "Negative")])

ggplot(data=threshold.results.b2.2$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b2.2$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))


## Repeat until all no negative cells remain (usually 3 rounds)...
## Round 3 /B1 
bar.table_sweep.list.b1 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b1[[n]] <- classifyCells(bar.table.b1, q=q)
  names(bar.table_sweep.list.b1)[n] <- paste("q=",q,sep="")
}

threshold.results.b1.3 <- findThresh(call.list=bar.table_sweep.list.b1)
b1.round3.calls <- classifyCells(bar.table.b1, q=findQ(threshold.results.b1.3$res, threshold.results.b1.3$extrema))
b1.neg.cells <- c(b1.neg.cells, names(b1.round3.calls)[which(b1.round3.calls == "Negative")])

ggplot(data=threshold.results.b1.3$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b1.3$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

## Repeat until all no negative cells remain (usually 3 rounds)...
b1.final.calls <- c(b1.round3.calls, rep("Negative",length(b1.neg.cells)))
names(b1.final.calls) <- c(names(b1.round3.calls),b1.neg.cells)
save(b1.final.calls, file="b1.final.calls.Rdata")


## Round 3 /B2
bar.table_sweep.list.b2 <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list.b2[[n]] <- classifyCells(bar.table.b2, q=q)
  names(bar.table_sweep.list.b2)[n] <- paste("q=",q,sep="")
}

threshold.results.b2.3 <- findThresh(call.list=bar.table_sweep.list.b2)
b2.round3.calls <- classifyCells(bar.table.b2, q=findQ(threshold.results.b2.3$res, threshold.results.b2.3$extrema))
b2.neg.cells <- c(b2.neg.cells, names(b2.round3.calls)[which(b2.round3.calls == "Negative")])

ggplot(data=threshold.results.b2.3$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() +
  geom_vline(xintercept=threshold.results.b2.3$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

## Repeat until all no negative cells remain (usually 3 rounds)...
b2.final.calls <- c(b2.round3.calls, rep("Negative",length(b2.neg.cells)))
names(b2.final.calls) <- c(names(b2.round3.calls),b2.neg.cells)
save(b2.final.calls, file="b2.final.calls.Rdata")


## Add back the prefix to cell IDs in B2 data to match cell ids in Seurat Object
names(b2.final.calls) <- paste0("2_", names(b2.final.calls))
final.calls <- c(b1.final.calls, b2.final.calls)


## Saving
exp@meta.data[,"MULTI"] <- "Unknown"
exp@meta.data[names(final.calls), "MULTI"] <- final.calls
saveRDS(exp, file = "exp.step1.rds")



exp@meta.data %>%group_by(MULTI) %>% summarise("number of cells" = n())
#Bar1: 7019 
#Bar2: 7512 
#Bar3: 7264 
#Bar4: 6895 
#Doublet: 4128 
#Negative: 2207 


