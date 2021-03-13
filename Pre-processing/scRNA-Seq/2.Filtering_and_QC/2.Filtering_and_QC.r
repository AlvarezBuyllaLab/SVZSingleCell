
# Redmond et al 2021 - Step 2: Filtering and QC"

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


#Removing doublets detected on demultiplex----
 exp <- readRDS("../1.Demultiplexing/exp.step1.rds")
#Bar1 and Bar2 are males/p35
#Bar 3 and 4, females/p29
Idents(exp) <- "MULTI"
exp <- subset(exp, idents = c("Bar1", "Bar2", "Bar3", "Bar4", "Negative"))


# Setting the cutoff for mitochondrial reads----
# We decided to work only with cells with less than 10% of reads mapping to the mitochondrial genome

exp[["percent.mt"]] <- PercentageFeatureSet(exp, pattern = "^mt-")

exp <- subset(exp, subset = percent.mt >0 &
                            percent.mt < 10)
#In the manuscript we used a low threshold of 0, which excluded 40 high-quality cells with 0 mitochondrial reads. Ideally, we would only set a high cutoff: 
#exp <- subset(exp, subset = percent.mt < 10)


#Setting the cutoffs for # of reads----

hist(exp@meta.data$nCount_RNA)
quantile(exp@meta.data$nCount_RNA)
countscutoff <- quantile(exp@meta.data$nCount_RNA, c(.05, .95))
exp <- subset(exp, subset = nCount_RNA > countscutoff[1] & 
                            nCount_RNA < countscutoff[2])


#Setting the cutoffs for # of genes----

hist(exp@meta.data$nFeature_RNA)
quantile(exp@meta.data$nFeature_RNA)
featscutoff <- quantile(exp@meta.data$nFeature_RNA, c(.05, .95))

exp <- subset(exp, subset = nFeature_RNA > featscutoff[1])

length(rownames(exp@meta.data)) #24261 cells


# Normalization, dimensionality reduction and clustering----
plan("multiprocess", workers = 24)
options(future.globals.maxSize= 600*1024^2)

exp <- SCTransform(exp, verbose = T)
exp <- RunPCA(exp, npcs=100, verbose = T)

exp <- RunUMAP(exp, reduction = "pca", dims = 1:100, verbose = T)

exp <- FindNeighbors(exp, dims = 1:100, verbose = T)
exp <- FindClusters(exp, resolution = c(0.8, seq(0.5, 2, 0.5)), verbose = T)

Idents(exp) <- "SCT_snn_res.1.5"
DimPlot(exp, label = T) + NoLegend() + coord_fixed()


# Saving
saveRDS(exp, file = "exp.step2.rds")
`
