exp <- readRDS("../2.Filtering_and_QC/exp.step2.rds")

options(future.globals.maxSize = 5000 * 1024^2)

LV.list <- SplitObject(exp, split.by = "well")
reference.list <- LV.list[c("well_1", "well_2")]
head(reference.list)
for (i in 1:length(x = LV.list)) {
    LV.list[[i]] <- SCTransform(object = LV.list[[i]], verbose = FALSE)
}

LV.features <- SelectIntegrationFeatures(object.list = LV.list)
LV.list <- PrepSCTIntegration(object.list = LV.list, anchor.features = LV.features, 
    verbose = FALSE)
LV.anchors <- FindIntegrationAnchors(object.list = LV.list, normalization.method = "SCT", 
    anchor.features = LV.features, verbose = FALSE)
exp.int <- IntegrateData(anchorset = LV.anchors, normalization.method = "SCT", 
    verbose = FALSE)

exp.int <- RunPCA(exp.int, verbose = T)
exp.int <- RunUMAP(exp.int, dims = 1:30, slot = "counts")
exp.int <- FindNeighbors(exp.int, dims = 1:30)
exp.int <-FindClusters(exp.int, resolution = c(0.8, seq(0.5, 2, 0.5)), verbose = T)

Idents(exp.int) <- "integrated_snn_res.1.5"
DefaultAssay(exp) <- "SCT"
saveRDS(exp, file = "exp.step3.rds")
