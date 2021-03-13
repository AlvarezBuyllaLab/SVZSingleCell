---
title: "Redmond et al 2021 - Demultiplexing"
---

#Setup

library(AUCell)
library(GSEABase)
library(Seurat)
library(shiny)
library(viridis)
library(scales)
library(scico)
library(ggplot2)
library(patchwork)
library(sctransform)
library(tidyr)


#Reading the Cellranger output and Creating a Seurat Object
```{r}
svz_matrix <- Read10X(
  data.dir = c("/data4/svz_10x/wholecell/data/pipestance/g1/outs/filtered_feature_bc_matrix",
               "/data4/svz_10x/wholecell/data/pipestance/g2/outs/filtered_feature_bc_matrix"),
  unique.features = TRUE,
  strip.suffix = FALSE
)

exp <- CreateSeuratObject(counts = svz_matrix$`Gene Expression`)
exp[['Protein']] <- CreateAssayObject(counts = svz_matrix$`Antibody Capture`)

exp

#Barcodes used: A3(TGAGACCT), A4(GCACACGC), A5(AGAGAGAG), A6 (TCACAGCA)
seurat.cells <- names(experiment$orig.ident)
cell.id.vec <- seurat.cells
multi.ref <- c("TGAGACCT", "GCACACGC", "AGAGAGAG", "TCACAGCA")



#1. Demultiplexing

##Loading FastQ files
```{r}
#tgacaat - B1
#acagtgat - B2

#B1:
readTable.b1 <- MULTIseq.preProcess(R1 = "/data4/svz_10x/wholecell/data/fastq/subsampled/TGACCAAT_S5_L001_R1_001_subsampled.fastq.gz", R2 = "/data4/svz_10x/wholecell/data/fastq/subsampled/TGACCAAT_S5_L001_R2_001_subsampled.fastq.gz", cellIDs = cell.id.vec)
save(readTable.b1, file= "readTable.b1.Rdata")

#B2:
readTable.b2 <- MULTIseq.preProcess(R1 = "/data4/svz_10x/wholecell/data/fastq/subsampled/ACAGTGAT_S6_L001_R1_001_subsampled.fastq.gz", R2 = "/data4/svz_10x/wholecell/data/fastq/subsampled/ACAGTGAT_S6_L001_R2_001_subsampled.fastq.gz", cellIDs = b2.cells.fixed)
save(readTable.b2, file= "readTable.b2.Rdata")
```

##Determining cell IDs in B1 and B2
```{r}
#This is a convoluted and non optimal way to determine the cell ids in B1 and B2. Instead, it should be determined using a filtering function to detect the b2 prefix "2_" as a regular expression and then separate the B1 and B2 ids. I am leaving this code here as this is how I originally determined B1 and B2 cells for future reference.

#all cells
all.cells <- names(experiment$orig.ident)

#B1:
b1.cells <- unique(readTable %>% select(Cell))
b1.cells <- as.vector(b1.cells[,1])
save(b1.cells, file = "b1.cells.Rdata")

#B2:
b2.cells <- all.cells[!(all.cells %in% b1.cells)]
save(b2.cells, file = "b2.cells.Rdata")

#b1.cells IDS:
b1.cells[1] #"ATTCCTACATACTGAC"
nchar(b1.cells[1]) #16 characters

b2.cells[1] #"2_AAACCCAAGAGCTGAC"
nchar(b2.cells[1]) #18 characters

#we have to remove the 2_ prefix of b2 ids in order to align reads:
b2.cells.fixed <- sub("2_", "", b2.cells)#b2.fixed was the reference used to align in the previous chunk 
```

##Checking the number of unique cells
```{r}
unique(readTable.b1 %>% select(Cell)) #18,028
unique(readTable.b2 %>% select(Cell)) #16,979
unique(readTable %>% select(Cell)) #18,046: the number of cells in G1. This is when we realized B2 wasn't being read.
```

## Perform MULTI-seq sample barcode alignment
```{r}
bar.table.b1 <- MULTIseq.align(readTable.b1, b1.cells, multi.ref)
save(bar.table.b1, file = "bar.table.b1.Rdata")

bar.table.b2 <- MULTIseq.align(readTable.b2, b2.cells.fixed, multi.ref)
save(bar.table.b2, file = "bar.table.b2.Rdata")
```

##Visual inspection in the barcode space
```{r}
#B1
bar.tsne.b1 <- barTSNE(bar.table.b1 [,1:4]) 
save(bar.tsne.b1, file="bar.tsne.b1.Rdata")

pdf("b1.bc.check.pdf")
for (i in 3:ncol(bar.tsne.b1)) {
    g <- ggplot(bar.tsne.b1, aes(x = TSNE1, y = TSNE2, color = bar.tsne.b1[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne.b1)[i]) +
    theme(legend.position = "none") 
    print(g)
}
dev.off()

#B2
bar.tsne.b2 <- barTSNE(bar.table.b2 [,1:4]) 
save(bar.tsne.b2, file="bar.tsne.b2.Rdata")

pdf("b2.bc.check.pdf")
for (i in 3:ncol(bar.tsne.b2)) {
    g <- ggplot(bar.tsne.b2, aes(x = TSNE1, y = TSNE2, color = bar.tsne.b2[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne.b2)[i]) +
    theme(legend.position = "none") 
    print(g)
}
dev.off()
```

## Round 1 / B1
```{r}
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
```

## Round 1 / B2 
```{r}
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
```

## Round 2 /B1
```{r}
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
```

## Round 2 /B2
```{r}
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
```

## Repeat until all no negative cells remain (usually 3 rounds)...
## Round 3 /B1 
```{r}
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
```

## Round 3 /B2
```{r}
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
```

## Add back the prefix to cell IDs in B2 data to incorporate in 
```{r}
names(b2.final.calls) <- paste0("2_", names(b2.final.calls))
final.calls <- c(b1.final.calls, b2.final.calls)

```

##Saving
```{r}
experiment.deMPx <- experiment
experiment.deMPx@meta.data[,"MULTI"] <- "Unknown"
experiment.deMPx@meta.data[names(final.calls), "MULTI"] <- final.calls
save(experiment.deMPx, file = paste0("WC_sctransformed_demultiplexed_",Sys.Date(),".Rdata"))
````