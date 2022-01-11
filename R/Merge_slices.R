library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
dir = c('data/Mouse_brain_A1-1/outs',
        'data/Mouse_brain_B1-1/outs',
        'data/Mouse_brain_C1-1/outs',
        'data/Mouse_brain_D1-1/outs',
		    'data/Mouse_brain_A1-2/outs',
        'data/Mouse_brain_B1-2/outs',
        'data/Mouse_brain_C1-2/outs',
        'data/Mouse_brain_D1-2/outs')
names(dir) = c('A1-1', 'B1-1', 'C1-1', 'D1-1','A1-2', 'B1-2', 'C1-2', 'D1-2')
brain <- list()
for(i in 1:length(dir)){
brain[[i]] <-Seurat::Load10X_Spatial(data.dir = dir[i])
brain[[i]]@meta.data$orig.ident <-names(dir)[i]
}
## SCTransform normalization
for (i in 1:length(brain)) {
    brain[[i]] <- SCTransform(brain[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)}
## merge合并
dir.create("Merge_Data")
setwd("/Merge_Data")
brain.merge <- merge(brain[[1]], y=c(brain[[2]], brain[[3]], brain[[4]],brain[[5]], brain[[6]], brain[[7]],brain[[8]]))
dim(brain.merge)
table(brain.merge@meta.data$orig.ident)
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain[[1]]), VariableFeatures(brain[[2]]), VariableFeatures(brain[[3]]), VariableFeatures(brain[[4]]),VariableFeatures(brain[[5]]), VariableFeatures(brain[[6]]), VariableFeatures(brain[[7]]), VariableFeatures(brain[[8]]))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge,resolution = 1.2, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
pdf("merge_umap.pdf",width=10,height=5)
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()
pdf("merge_umap_slice.pdf",width=40,height=5)
SpatialDimPlot(brain.merge)
dev.off()
brain.merge@meta.data$new_clusters <-paste(brain.merge$seurat_clusters,brain.merge$orig.ident,sep = "_")
saveRDS(brain.merge,"merge.rds") # Save the workspace

## Find marker genes, setting logfc.threshold = 0.25
brain.markers <- FindAllMarkers(brain.merge, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
write.table(brain.markers,"marker.txt",row.names=TRUE,col.names=TRUE,sep="\t")
write.table(brain.merge@meta.data,"meta.data.txt",row.names=TRUE,col.names=TRUE,sep="\t")

topgene<-brain.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf("merge_heatmap.pdf",width=12,height=8)
DoHeatmap(brain.merge, features = topgene$gene,size = 2) + NoLegend()
dev.off()

#gene expression
dir.create("Gene_exp")
setwd("Gene_exp")
#marker genes top30 Vio
for(i in 0:(length(unique(brain.markers[,"cluster"]))-1)){
ind <-which(brain.markers[,"cluster"]==i)[1:30]
gene <-brain.markers[ind,"gene"]
pdf(paste("Cluster",i,"Vio.pdf",sep ="_"),width=ceiling(length(unique(brain.markers[,"cluster"]))/3),height=length(unique(brain.markers[,"cluster"]))*3)
print(VlnPlot(brain.merge,features =gene, group.by = "seurat_clusters",pt.size = 0,ncol = 1))
dev.off()
}
#marker genes top10 Umap
for(i in 0:(length(unique(brain.markers[,"cluster"]))-1)){
ind <-which(brain.markers[,"cluster"]==i)[1:10]
gene <-brain.markers[ind,"gene"]
pdf(paste("Cluster",i,"Umap.pdf",sep ="_"),width=10,height=25)
print(FeaturePlot(brain.merge, features =gene, cols = c("grey", "red"),reduction = "umap",ncol =2))
dev.off()
}
####marker genes top10 Spatial
for(i in 0:(length(unique(brain.markers[,"cluster"]))-1)){
ind <-which(brain.markers[,"cluster"]==i)[1:10]
gene <-brain.markers[ind,"gene"]
pdf(paste("Cluster",i,"Spatial.pdf",sep ="_"),width=10,height=30)
print(SpatialFeaturePlot(brain.merge, features = gene,ncol =1 ))
dev.off()
