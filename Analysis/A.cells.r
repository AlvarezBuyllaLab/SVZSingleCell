---
#title: "Cebrian-Silla, Nascimento, Redmond, Mansky et al. 2021 - Analysis: A Cells"

---
#setup

library(knitr)
library(future)
library(ggplot2) 
library(reticulate) 
library(Seurat)
library(R.utils)
library(dplyr)
library(tibble)
library(viridis)
library(patchwork)
library(dendextend)
library(stats)
library(remotes)
library(readr)
library(biomaRt)
library(enrichplot)
library(scales)
library(tidyr)
library(scico)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Loading the annotated Seurat Object

exp <- readRDS("../Files/S")
Idents(exp) <- "integrated_snn_res.1.5"
DimPlot(exp, label = T) + coord_fixed() + NoLegend()

exp.lin <- readRDS("/data4/marcos/WC_Mouse.2020/exp.lin.rds")
DimPlot(exp.lin, label = T) + coord_fixed() + NoLegend()

##Testing a single gene expression:

FeaturePlot(exp, "Vax1", order = T) + scale_color_viridis(option = "A") + coord_fixed()


#Find All Markers
plan("multiprocess", workers = 16)
plan()
all.markers.1.5 <- FindAllMarkers(exp)
saveRDS(all.markers.1.5, file = "all.markers.1.5.rds")
#A Cells clusters are: 6, 4, 1 and 0


##Listing markers for each cluster
cluster15.allmarkers <- all.markers.1.5 %>% filter(cluster==15 & avg_logFC >0) %>% pull(gene)
cluster6.allmarkers <- all.markers.1.5 %>% filter(cluster==6 & avg_logFC >0) %>% pull(gene)
cluster4.allmarkers <- all.markers.1.5 %>% filter(cluster==4 & avg_logFC >0) %>% pull(gene)
cluster1.allmarkers <- all.markers.1.5 %>% filter(cluster==1 & avg_logFC >0) %>% pull(gene)
cluster0.allmarkers <- all.markers.1.5 %>% filter(cluster==0 & avg_logFC >0) %>% pull(gene)

write.table(cluster15.allmarkers, quote = F, row.names = F, col.names = F, "cluster15.allmarkers.txt", sep=",")
write.table(cluster6.allmarkers, quote = F, row.names = F, col.names = F, "cluster6.allmarkers.txt", sep=",")
write.table(cluster4.allmarkers, quote = F, row.names = F, col.names = F, "cluster4.allmarkers.txt", sep=",")
write.table(cluster1.allmarkers, quote = F, row.names = F, col.names = F, "cluster1.allmarkers.txt", sep=",")
write.table(cluster0.allmarkers, quote = F, row.names = F, col.names = F, "cluster0.allmarkers.txt", sep=",")


##Plots using all markers 
#cluster6
p1 <- FeaturePlot(exp.lin, features = cluster6.allmarkers[1:9], combine = FALSE, pt.size = 0.2, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() +  fix.sc + simple + xlim(-6, 7) + ylim(-15, 8))

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 6 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster6.all.markers.png", width = 15, height = 15, units = "cm", scale = 2)


#cluster4
p1 <- FeaturePlot(exp, features = cluster4.allmarkers[1:9], combine = FALSE, pt.size = 0.2, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 4 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster4.all.markers.png", width = 15, height = 15, units = "cm", scale = 2)


#cluster1
p1 <- FeaturePlot(exp, features = cluster1.allmarkers[1:9], combine = FALSE, pt.size = 0.2, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 1 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster1.all.markers.png", width = 15, height = 15, units = "cm", scale = 2)


#cluster0
p1 <- FeaturePlot(exp, features = cluster0.allmarkers[1:9], combine = FALSE, pt.size = 0.2, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 0 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster0.all.markers.png", width = 15, height = 15, units = "cm", scale = 2)
```
##GO analysis
```{r}
#visago bioconductor package to do GO
#clusterprofiler another package to plot GO
```


```{r}
#GO analysis was performed by importing the markers from each cluster into pantherdb, using the .txt files created in chunk#5 of this notebook.

#A for loop to read all .txt files exported from pantherdb and create tables names "cluster#.go"
clusters <- c(0,1,4,6,15)
for (i in clusters) {
  assign(paste0("cluster",i,".go"), read_delim(paste0("GO/cluster", i,  ".allmarkers.txt"), 
    "\t", escape_double = FALSE, col_names = c("go", "ref", "obs", "expected", "pos.neg", "fold.enrich", "raw.p.value", "fdr"), 
    trim_ws = TRUE, skip = 12))
} 
clusterorder <- c("Cluster 1", "Cluster 0", "Cluster 4", "Cluster 6", "Cluster 15")

#create a dataframe with all categories and clusters:
go.cat.n <- nrow(cluster1.go)
go.data <- do.call("rbind", list(cluster1.go, cluster0.go, cluster4.go, cluster6.go, cluster15.go))
go.data[1:go.cat.n,"cluster"] <- "Cluster 1"
go.data[(1+go.cat.n):(2*go.cat.n),"cluster"] <- "Cluster 0"
go.data[(1+(2*go.cat.n)):(3*go.cat.n),"cluster"] <- "Cluster 4"
go.data[(1+(3*go.cat.n)):(4*go.cat.n),"cluster"] <- "Cluster 6"
go.data[(1+(4*go.cat.n)):(5*go.cat.n),"cluster"] <- "Cluster 15"
go.data$cluster <- factor(go.data$cluster, levels = clusterorder)

saveRDS(go.data, file = "A.cells.go.data.rds")
                   
#plotting barplot
df <- go.data %>% filter(raw.p.value < 0.01 & fdr <0.05 & cluster =="Cluster 0") %>% top_n(20, fold.enrich)


df$go <- factor(df$go, levels= rev(unique(as.character(df$go))))
#df <- transform(df, go=reorder(variable, -value)) #to reorder the rows according to another value

#single cluster
ggplot(df, aes(x= fold.enrich, y = go, fill = fdr)) + 
  scale_fill_gradient() + 
  geom_bar(stat = "identity") + 
  scale_fill_continuous(high = "#fff7ec", low = "#7f0000") + 
  theme_minimal() + 
  theme(panel.grid.major.y = element_blank()) + 
  labs(title = "Cluster 0")

#allclusters wrapped
##Finding the top N categories enriched in each cluster, filtering for p value < 0.01 and FDR <0.05:
go.cats <- unique(go.data %>% filter(raw.p.value < 0.01 & fdr <0.05) %>% group_by(cluster) %>% slice_max(n = 10, order_by = fold.enrich, with_ties = F) %>% pull(go))
go.cats.curated <- c("mitotic DNA replication initiation (GO:1902975)",
                     "DNA strand elongation involved in DNA replication (GO:0006271)",
                     "DNA replication, synthesis of RNA primer (GO:0006269)",
                     "mitotic DNA replication (GO:1902969)",
                     "positive regulation of DNA-directed DNA polymerase activity (GO: 1900264)",
                     "regulation of postsynaptic density protein 95 clustering (GO:1902897)",
                      "nucleosome disassembly (GO:0006337)",
                     "ephrin receptor signaling pathway (GO:0048013)",
                     "commitment of multipotent stem cells to neuronal lineage in forebrain (GO:0021898)",
                     "axon guidance (GO:0007411)",
                     "chemotaxis (GO:0006935)",
                    "regulation of postsynaptic density organization (GO:1905874)",
                     "neural precursor cell proliferation (GO:0061351)",
                     "olfactory bulb interneuron differentiation (GO:0021889)",
                     "regulation of synaptic vesicle clustering (GO:2000807)",
                     "regulation of negative chemotaxis (GO:0050923)",
                     "neuroligin clustering involved in postsynaptic membrane assembly (GO:0097118)",
                     "positive regulation of dendritic spine morphogenesis (GO:0061003)",
                     "cerebral cortex radial glia guided migration (GO:0021801)",
                     "anterograde dendritic transport (GO:0098937)",
                     "modulation of inhibitory postsynaptic potential (GO:0098828)",
                     "spontaneous synaptic transmission (GO:0098814)",
                     "spontaneous neurotransmitter secretion (GO:0061669)")

go.cats.curated <- rev(unique(go.data %>% filter(go %in% go.cats.curated) %>% group_by(cluster) %>% top_n(30, -fdr) %>% pull(go)))

#creating dataframe to use in plots:
df.all<- go.data %>% filter(go %in% go.cats.curated)
df.all <- df.all %>% mutate(fold.enrich = replace_na(fold.enrich, 0))
df.all <- df.all %>% mutate(sig = ifelse(fdr < 0.05, T, F))
df.all$cluster <- factor(df.all$cluster, levels = clusterorder)
df.all$go <- factor(df.all$go, levels = go.cats.curated)

#Heatmap
ggplot(df.all %>% mutate(fold.enrich = replace_na(fold.enrich, 0)), aes(cluster, go)) + geom_tile(aes(fill = fold.enrich)) + scale_fill_viridis() + theme_minimal()

#Dotplot
ggplot(df.all, aes(cluster, go)) + 
  geom_point(aes(col = sig, size = fold.enrich)) + 
  scale_color_manual(name = NULL, labels = c("FDR p > 0.05", "FDR p < 0.05"), values = viridis(2, direction = -1)) + 
  scale_size(name = "Fold enrichment") +
  scale_x_discrete(labels=c("Cluster 1" = "1", "Cluster 0" = "0", "Cluster 4" = "4", "Cluster 6" = "6", "Cluster 15" = "15")) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), axis.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal")

ggsave(filename = "go.dotplot.png", width = 3.6, height = 2.1, scale =2)

```

#Subsetting A cells
```{r}
exp.a <-subset(exp, idents = c(15,6,1,4,0))

a.cells.palette <- c("#671313", "#9B1C1C", "#D93030", "#DF5353", "#EC9898")

a.cells <- WhichCells(exp, idents = c(15,6,1,4,0))
exp@meta.data[a.cells,"acells"] <- exp@meta.data[a.cells,"integrated_snn_res.1.5"]


DimPlot(exp, group.by = "acells", label=T, shuffle = T, cols = a.cells.palette ) + 
  coord_fixed() + 
  xlim(-6, 7) + ylim(-15, 8) + 
  NoAxes() + 
  labs(title = "A Cells")
ggsave

DimPlot(exp, cells = a.cells, label=T, shuffle = T, cols = a.cells.palette ) + 
  xlim(-4.5, 3.5) + ylim(-14, -2) + 
  coord_fixed() + 
  NoAxes() + 
  labs(title = "A Cells Subset")

DimPlot(exp, cells = a.cells, label=F, group.by = "Phase") + 
  xlim(-4.5, 3.5) + ylim(-14, -2) + 
  coord_fixed() + 
  NoAxes() + 
  labs(title = "A Cells - Cell cycle labeling")

a.markers <- FindAllMarkers(exp.a, only.pos = T)
saveRDS(exp.a, file="exp.a.RDS")
saveRDS(a.markers, file="a.markers.RDS")
```

##Checking gene expression
###Avg Expression
```{r}
Idents(exp) <- "integrated_snn_res.1.5"
exp.avg <- AverageExpression(exp)

exp.avg$integrated %>% top_n(20)
```

###Plotting
```{r}
gene = "Slit1"

a <- FeaturePlot(exp.a, features = gene, order = T, pt.size = 0.5) + 
  scale_color_viridis(alpha = 1, option = "A") + 
  coord_fixed() +
  xlim(-4.5, 3.5) + 
  ylim(-14, -2) +
  simple

b <- VlnPlot(exp.a, features = gene, pt.size = 0) + NoLegend()
a+b

VlnPlot(exp.a, features = c(""))
```
#exp.a Markers
##FindAllMarkers in exp.a
```{r}
cutoff <- 0.2  #using 0.2 leads to a better defined, shorter list

cluster15.markers <- a.markers %>% filter(cluster==15 & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
cluster6.markers <- a.markers %>% filter(cluster==6 & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
cluster4.markers <- a.markers %>% filter(cluster==4 & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
cluster1.markers <- a.markers %>% filter(cluster==1 & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
cluster0.markers <- a.markers %>% filter(cluster==0 & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)

write.table(cluster15.markers, quote = F, row.names = F, col.names = F, "cluster15.strictmarkers.txt", sep=",")
write.table(cluster6.markers, quote = F, row.names = F, col.names = F, "cluster6.strictmarkers.txt", sep=",")
write.table(cluster4.markers, quote = F, row.names = F, col.names = F, "cluster4.strictmarkers.txt", sep=",")
write.table(cluster1.markers, quote = F, row.names = F, col.names = F, "cluster1.strictmarkers.txt", sep=",")
write.table(cluster0.markers, quote = F, row.names = F, col.names = F, "cluster0.strictmarkers.txt", sep=",")

#Judging by the number of markers, clusters 0 (dorsal), 1 (ventral), and 6 (proliferative) are very unique. Cluster 4 seems like a transitional identity

#String Analysis:
  #Cluster1/Ventral: 
    # 0.2 cutoff: https://version-11-0.string-db.org/cgi/network.pl?networkId=WWV16gHY9tKm #20 genes
    # no cutoff: https://version-11-0.string-db.org/cgi/network.pl?networkId=rTkKFUYrPbYl #121 genes

  #Cluster0/Dorsal: 
    #https://version-11-0.string-db.org/cgi/network.pl?networkId=8WC9g7rgGcxQ #10 genes
    #https://version-11-0.string-db.org/cgi/network.pl?networkId=OHEYxcivrclk #60 genes

  #Cluster6: Early neuroblasts

"Robo1" %in% cluster1.markers
```

####Plots
```{r}
#cluster 6
p1 <- FeaturePlot(exp.a, features = cluster6.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) + fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 6 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster6.markers.png", width = 15, height = 15, units = "cm", scale = 2)

#cluster 4
p1 <- FeaturePlot(exp.a, features = cluster4.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 4 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("cluster4.markers.png", width = 15, height = 15, units = "cm", scale = 2)


#cluster 1
p1 <- FeaturePlot(exp.a, features = cluster1.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 1 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))
ggsave("cluster1.markers.png", width = 15, height = 15, units = "cm", scale = 2)


#cluster 0
p1 <- FeaturePlot(exp.a, features = cluster0.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 0 markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))
ggsave("cluster0.markers.png", width = 15, height = 15, units = "cm", scale = 2)

```
###Cluster 1 vs Cluster 0
```{r}
Idents(exp.a) <- "integrated_snn_res.1.5"
vxd.a.markers <- FindMarkers(exp.a, ident.1 = 1, ident.2 = 0)

v.a.markers <- rownames(vxd.a.markers %>% filter(avg_logFC >0 & pct.2 < 0.5))
d.a.markers <- rownames(vxd.a.markers %>% filter(avg_logFC <0 & pct.1 < 0.5))

write.table(d.a.markers, quote = F, row.names = F, col.names = F, "d.a.markers.txt", sep=",")
```

####Plots
```{r}
#cluster 1
Idents(exp.a) <- "cellid"

##Feat.
p1 <- FeaturePlot(exp.a, features = v.a.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 1 markers",
    subtitle = "Identified comparing clusters 1 and 0",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("ventral.a.markers.feat.png", width = 15, height = 15, units = "cm", scale = 2)

VlnPlot(exp.a, features = v.a.markers[1:20], stack = T, flip = T, idents = c("Ventral","Dorsal")) + NoLegend()

##Stacked Violins
StackedVlnPlot(obj = exp.a, features = v.a.markers[1:20], idents = c("Ventral","Dorsal"))
ggsave("ventral.a.markers.vln.png", width = 7, height = 15, units = "cm", scale = 2)
StackedVlnPlot(obj = exp.a, features = v.a.markers[21:40], idents = c("Ventral","Dorsal"))
ggsave("ventral.a.markers.vln2.png", width = 7, height = 15, units = "cm", scale = 2)

#cluster 0
p1 <- FeaturePlot(exp.a, features = d.a.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top cluster 0 markers",
    subtitle = "Identified comparing clusters 1 and 0",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))
ggsave("dorsal.a.markers.feat.png", width = 15, height = 15, units = "cm", scale = 2)

```



#TestPlots
```{r}
FeaturePlot(exp.a, features = "Nxph1", pt.size = 1, order = T) + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) +  fix.sc + simple

DimPlot(exp.a, group.by = "integrated_snn_res.0.8", label = T) + coord_fixed()
DimPlot(exp.a, group.by = "integrated_snn_res.1.5", label = T) + coord_fixed()

"Vax1" %in% v.a.markers
```


#Assigning A cell identities
```{r}
Idents(exp.a) <- "integrated_snn_res.1.5"
v.acells <- WhichCells(exp.a, idents = 1)
d.acells <- WhichCells(exp.a, idents = 0)
nb.acells <- WhichCells(exp.a, idents = c(6, 15))
inter.acells <- WhichCells(exp.a, idents = 4)


exp.a@meta.data[c(v.acells, inter.acells, d.acells),  "cellid"] <- "Young Neurons"
exp.a@meta.data[nb.acells,  "cellid"] <- "Neuroblasts"


Idents(exp.a) <- "cellid"
#Idents(exp.a) <- factor(Idents(exp.a),levels = c("B ventral", "A ventral", "B dorsal", "A dorsal"))

DimPlot(exp.a, group.by = "cellid", label = T) + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -2) + simple

a.id.markers <- FindAllMarkers(exp.a, only.pos = T)
```

##Markers
```{r}
cutoff <- 1  #using 0.2 leads to a better defined, shorter list

ventral.id.markers <- a.id.markers %>% filter(cluster== "Ventral" & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
dorsal.id.markers <- a.id.markers %>% filter(cluster== "Dorsal" & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
nb.id.markers <- a.id.markers %>% filter(cluster=="Neuroblast" & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)
inter.id.markers <- a.id.markers %>% filter(cluster=="Transition" & avg_logFC >0 & pct.2 < cutoff) %>% pull(gene)


write.table(ventral.id.markers, quote = F, row.names = F, col.names = F, "ventral.id.markers.txt", sep=",")
write.table(dorsal.id.markers, quote = F, row.names = F, col.names = F, "dorsal.id.markers.txt", sep=",")
write.table(nb.id.markers, quote = F, row.names = F, col.names = F, "nb.id.markers.txt", sep=",")
write.table(inter.id.markers, quote = F, row.names = F, col.names = F, "inter.id.markers.txt", sep=",")

#Judging by the number of markers, clusters 0 (dorsal), 1 (ventral), and 6 (proliferative) are very unique. Cluster 4 seems like a transitional identity
```

##Feature Plots
```{r}
#Ventral
p1 <- FeaturePlot(exp.a, features = ventral.id.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) + fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top ventral young neuron markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("ventral.id.markers.png", width = 15, height = 15, units = "cm", scale = 2)

#Dorsal
p1 <- FeaturePlot(exp.a, features = dorsal.id.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) + fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top dorsal young neuron markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("dorsal.id.markers.png", width = 15, height = 15, units = "cm", scale = 2)

#Neuroblast
p1 <- FeaturePlot(exp.a, features = early.id.markers[1:9], combine = FALSE, pt.size = 0.6, order = T)
fix.sc <- scale_color_viridis(alpha = 1, option = "A") 
simple <- NoAxes() + NoLegend() 
Plot1 <- lapply(p1, function (x) x + coord_fixed() + xlim(-4.5, 3.5) + ylim(-14, -5) + fix.sc + simple)

wrap_plots(Plot1, nrow = 3) + 
  plot_annotation(
    title = "Top  neuroblast markers",
    subtitle = "Identified using FindAllMarkers()",
    caption = paste0("Generated on ", format(Sys.time(), "%D"), ". Code available on A.cells.Rmd"))

ggsave("early.id.markers.png", width = 15, height = 15, units = "cm", scale = 2)
```
##Violin Plots
```{r}
#Dorsal
VlnPlot(exp.a, features = dorsal.id.markers[1:9], pt.size = 0) + NoLegend()

#Ventral
VlnPlot(exp.a, features = ventral.id.markers[1:9], pt.size = 0) + NoLegend()

#Early
VlnPlot(exp.a, features = early.id.markers[1:9], pt.size = 0) + NoLegend()
```
#Heatmap + Dendrogram
```{r}
exp.a.heatmap <- SCTransform(exp.a, new.assay.name = "SCT2", variable.features.n = 20681, vars.to.regress = "percent.mt")
exp.a.heatmap <- BuildClusterTree(exp.a.heatmap, reorder=T)

#Dendrogram
dend1 <- Tool(exp.a.heatmap, 'BuildClusterTree')
dend1
dend1 <- as.dendrogram(dend1)
ggplot(dend1)
ggsave("acells.dendrogram.png", scale =2, width = 4, height = 0.75)

mge.markers <- c("Gad1", "Gad2", "Lhx6", "Mapt", "Mllt11", "Myt1l", "Sez612")
lge.markers <- c("Meis2", "Ebf1", "Pcp4", "Isl1")

#Heatmap
heatmap<- DoHeatmap(subset(exp.a.heatmap, downsample = 200), features = rev(c(cluster1.markers[1:5], cluster0.markers[1:5], cluster4.markers[1:5], cluster6.markers[1:5],cluster15.markers[1:5])), size = 3)
ggsave("acells.heatmap.0.2threshold.png", scale =2, width = 4, height = 3)
```

#Feature plots for reginal markers
```{r}
FeaturePlot(exp.a, features = "Slit2", order = T, pt.size = 0.7) + 
  scale_color_viridis(alpha = 1, option = "A") + 
  coord_fixed() +
  xlim(-4.5, 3.5) + 
  ylim(-14, -2) +
  simple
ggsave(filename = "slit2.png", height = 2, width = 2, scale = 2)

FeaturePlot(exp.a, features = "Vax1", order = T, pt.size = 0.7) + 
  scale_color_viridis(alpha = 1, option = "A") + 
  coord_fixed() +
  xlim(-4.5, 3.5) + 
  ylim(-14, -2) +
  simple
ggsave(filename = "vax1.png", height = 2, width = 2, scale = 2)

FeaturePlot(exp.a, features = "Pax6", order = T, pt.size = 0.7) + 
  scale_color_viridis(alpha = 1, option = "A") + 
  coord_fixed() +
  xlim(-4.5, 3.5) + 
  ylim(-14, -2) +
  simple
ggsave(filename = "pax6.png", height = 2, width = 2, scale = 2)

FeaturePlot(exp.a, features = "Rlbp1", order = T, pt.size = 0.7) + 
  scale_color_viridis(alpha = 1, option = "A") + 
  coord_fixed() +
  xlim(-4.5, 3.5) + 
  ylim(-14, -2) +
  NoAxes()
ggsave(filename = "rlbp1.png", height = 2, width = 2, scale = 2)

```


#Trying enrichplot/clusterprofiler
#Vignette
```{r}
data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="ALL", pAdjustMethod = "fdr", readable=TRUE)

goplot(ego) #doesnt work for some reason

barplot(ego, showCategory=20)
dotplot(ego, showCategory=30)
dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

```
#ACell
```{r}
cluster.15.entrez <- bitr(cluster15.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
ego.15 <- enrichGO(cluster.15.entrez, OrgDb = "org.Mm.eg.db", ont="BP", pAdjustMethod = "fdr", readable=TRUE)

cluster.6.entrez <- bitr(cluster6.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
ego.6 <- enrichGO(cluster.6.entrez, OrgDb = "org.Mm.eg.db", ont="BP", pAdjustMethod = "fdr", readable=TRUE)

cluster.4.entrez <- bitr(cluster4.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
ego.4 <- enrichGO(cluster.4.entrez, OrgDb = "org.Mm.eg.db", ont="BP", pAdjustMethod = "fdr", readable=TRUE)

cluster.1.entrez <- bitr(cluster1.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
ego.1 <- enrichGO(cluster.1.entrez, OrgDb = "org.Mm.eg.db", ont="BP", pAdjustMethod = "fdr", readable=TRUE)

cluster.0.entrez <- bitr(cluster0.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
ego.0 <- enrichGO(cluster.0.entrez, OrgDb = "org.Mm.eg.db", ont="BP", pAdjustMethod = "fdr", readable=TRUE)
```

```{r}
dotplot(ego.15, showCategory=30) + labs(title = "Cluster 15")
dotplot(ego.6, showCategory=30) + labs(title = "Cluster 6")
dotplot(ego.4, showCategory=30) + labs(title = "Cluster 4")
dotplot(ego.1, showCategory=30) + labs(title = "Cluster 1")
dotplot(ego.0, showCategory=30) + labs(title = "Cluster 0")
```

