
library(tidyverse)
library(Seurat)
#library(fgsea)
library(cowplot)

#import male LHA data. n = 3 males, but pooled by authors into one sample run

male.lha.data <- Read10X(data.dir = "~/Box/Lateral Hypothalamus/LHA/male/")
female.lha.data <- Read10X(data.dir = "~/Box/Lateral Hypothalamus/LHA/female/")

male <- CreateSeuratObject(counts = male.lha.data, project = "maleLHA", min.cells = 3, min.features = 200)
female <- CreateSeuratObject(counts = female.lha.data, project = "femaleLHA", min.cells = 3, min.features = 200)

male <- PercentageFeatureSet(object = male, pattern = "^mt-", col.name = "percent.mt")
female <- PercentageFeatureSet(object = female, pattern = "^mt-", col.name = "percent.mt")

VlnPlot(object = male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
VlnPlot(object = female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

male <- subset(x = male, subset = nFeature_RNA > 1000 & percent.mt < 15)
female <- subset(x = female, subset = nFeature_RNA > 1000 & percent.mt < 15)


table(male@meta.data$orig.ident)

male <- male %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

female <- female %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

lha.anchors <- FindIntegrationAnchors(object.list = list(male, female), dims = 1:20)
lha.combined <- IntegrateData(anchorset = lha.anchors, dims = 1:20)

table(lha.combined@meta.data$orig.ident)


DefaultAssay(lha.combined) <- "integrated"

# Clustering. I have not toyed with these paramaters much, just adjusted res lower to grab neurons easier
lha.combined <- ScaleData(lha.combined, verbose = FALSE)
lha.combined <- RunPCA(lha.combined, npcs = 30, verbose = FALSE)

ElbowPlot(lha.combined)

lha.combined <- RunUMAP(lha.combined, reduction = "pca", dims = 1:10)

lha.combined <- FindNeighbors(lha.combined, reduction = "pca", dims = 1:10)
lha.combined <- FindClusters(lha.combined, resolution = 0.2)

#Neuronal markers. Looks like clusters 0, 1, 6, 10, and maybe 5 (probably not - only by meg3) are neurons.
#actually most neuronal clusters aren't showing glut or gaba markers. plus bad clusters
FeaturePlot(lha.combined, reduction = "umap", cols = c("lightgrey","black"), features = c("Snap25", "Syp", "Tubb3", "Elavl2", "Meg3"), pt.size = .1)
p2 <- DimPlot(lha.combined, reduction = "umap", label = TRUE) + NoLegend()

plot_grid(p1, p2)

DimPlot(lha.combined, reduction = "umap", group.by = "orig.ident")
#gaba / glut
FeaturePlot(lha.combined, reduction = "umap", cols = c("lightgrey","black"), features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"), pt.size = .1)



#extract neurons
lha.neurons <- subset(lha.combined, ident = c("0", "1", "6", "10"))

#re-cluster

lha.neurons <- SCTransform(lha.combined, variable.features.n = 3000)

lha.neurons <- RunPCA(lha.neurons, npcs = 30, verbose = FALSE)

VariableFeaturePlot(lha.neurons)
ElbowPlot(lha.neurons)

lha.neurons <- RunUMAP(lha.neurons, reduction = "pca", dims = 1:15)

lha.neurons <- FindNeighbors(lha.neurons, reduction = "pca", dims = 1:15)
lha.neurons <- FindClusters(lha.neurons, resolution = 0.2)

DimPlot(lha.neurons, reduction = "umap", label = TRUE) + NoLegend()

#gaba / glut
FeaturePlot(lha.neurons, reduction = "umap", cols = c("lightgrey","black"), features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"), pt.size = .1)

#sst in all neurons
FeaturePlot(lha.neurons, reduction = "umap", cols = c("lightgrey","black"), features = c("Sst"), pt.size = .1)

#needs more filtering I believe
lha.neurons.sub <- subset(lha.neurons, ident = c("0", "1", "3", "13"))

lha.neurons.sub <- SCTransform(lha.neurons.sub, variable.features.n = 3000)

lha.neurons.sub <- RunPCA(lha.neurons.sub, npcs = 30, verbose = FALSE)

VariableFeaturePlot(lha.neurons.sub)
ElbowPlot(lha.neurons.sub)

lha.neurons.sub <- RunUMAP(lha.neurons.sub, reduction = "pca", dims = 1:18)

lha.neurons.sub <- FindNeighbors(lha.neurons.sub, reduction = "pca", dims = 1:18)
lha.neurons.sub <- FindClusters(lha.neurons.sub, resolution = 1.6)

DimPlot(lha.neurons.sub, reduction = "umap", label = TRUE) + NoLegend()
plot_grid(p1, p3)
#gaba / glut
FeaturePlot(lha.neurons.sub, reduction = "umap", cols = c("lightgrey","black"), features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"), pt.size = .1)
#looks pretty good now

p2 <- FeaturePlot(lha.neurons.sub, reduction = "umap", cols = c("lightgrey","black"), features = c("Sst"), pt.size = .1)

plot_grid(p1, p2)

table(lha.neurons.sub@meta.data$orig.ident)

neuron.markers <- FindAllMarkers(lha.neurons.sub, only.pos = T)

top_markers(neuron.markers)
testing <- rename_clusters(lha.neurons.sub, marker.list = neuron.markers)

DimPlot(testing, reduction = "umap", label = TRUE) + NoLegend()


#sst clusters are 10, 17, 21
lha.cluster.degs <- Cluster_DEGs(lha.neurons.sub, condition.1 = "femaleLHA", condition.2 = "maleLHA")
c10 <- lha.cluster.degs[[11]]
c17 <- lha.cluster.degs[[18]]
c21 <- lha.cluster.degs[[22]]

BellagioPlot(lha.cluster.degs)

neuron.markers %>%
  mutate(sig = ifelse(p_val_adj < .05, "Y", "N")) %>%
  arrange(desc(sig), desc(avg_logFC)) %>%
  head(10)

test <- top_markers(neuron.markers, 10)


test <- neuron.markers %>%
  mutate(pi = avg_logFC * -log10(p_val))

#decrease res to zero to see all neurons
lha.neurons.sub <- FindClusters(lha.neurons.sub, resolution = 0)
DimPlot(lha.neurons.sub, reduction = "umap", label = TRUE) + NoLegend()

lha.allneuron.degs <- Cluster_DEGs(lha.neurons.sub, condition.1 = "femaleLHA", condition.2 = "maleLHA")

t <- lha.allneuron.degs[[1]]



#can we easily get degs for all glut or gaba? Or do we have to extract? looks like would have to extract
lha.neurons.sub <- FindClusters(lha.neurons.sub, resolution = 0.2)
DimPlot(lha.neurons.sub, reduction = "umap", label = TRUE) + NoLegend()
FeaturePlot(lha.neurons.sub, reduction = "umap", cols = c("lightgrey","black"), features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"), pt.size = .1)


Cluster_DEGs <- function(object, condition.1, condition.2, logfc.threshold = .1, min.pct = .05, test.use = "wilcox", logFCcollapse = 100){
  tryCatch({
    num.clusters <- nlevels(object@active.ident)
    Idents(object) <- paste0(object@active.ident, "_", object@meta.data$orig.ident)
    temp <- list()
    for(i in 1:num.clusters){
      tryCatch({
        x <- paste0(i-1, "_", condition.1)
        y <- paste0(i-1, "_", condition.2)
        temp[[i]] <- FindMarkers(object, ident.1 = x, ident.2 = y, logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use)
        temp[[i]] <- rownames_to_column(temp[[i]], var = "gene")
        temp[[i]]$avg_logFC[temp[[i]]$avg_logFC < -logFCcollapse] <- -logFCcollapse
        temp[[i]]$avg_logFC[temp[[i]]$avg_logFC > logFCcollapse] <- logFCcollapse
        temp[[i]][,"Cluster"] <- (i-1)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    return(temp)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
