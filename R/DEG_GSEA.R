library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
brain <-readRDS("/Merge_data/merge.rds")
#differential expression analysis (DEA) 
#C2 Fiber
dir.create("C2")
setwd("C2")
DefaultAssay(brain) <- "Spatial"
brain <- NormalizeData(brain)
brain@active.ident <-factor(brain$new_clusters)
marker_A <-FindMarkers(brain,ident.1="2_A1-2",ident.2="2_A1-1",logfc.threshold=log2(1.1),assay = "Spatial",slot = 'data')
marker_A <-marker_A[which(marker_A[,1]<=0.05),]
up_A <-rownames(marker_A[which(marker_A[,2]>=0),])
down_A <-rownames(marker_A[which(marker_A[,2]<=0),])
length(up_A)
length(down_A)
marker_B <-FindMarkers(brain,ident.1="2_B1-2",ident.2="2_B1-1",logfc.threshold=log2(1.1),assay = "Spatial",slot = 'data')
marker_B <-marker_B[which(marker_B[,1]<=0.05),]
up_B <-rownames(marker_B[which(marker_B[,2]>=0),])
down_B <-rownames(marker_B[which(marker_B[,2]<=0),])
length(up_B)
length(down_B)
marker_C <-FindMarkers(brain,ident.1="2_C1-2",ident.2="2_C1-1",logfc.threshold=log2(1.1),assay = "Spatial",slot = 'data')
marker_C <-marker_C[which(marker_C[,1]<=0.05),]
up_C <-rownames(marker_C[which(marker_C[,2]>=0),])
down_C <-rownames(marker_C[which(marker_C[,2]<=0),])
length(up_C)
length(down_C)
marker_D <-FindMarkers(brain,ident.1="2_D1-2",ident.2="2_D1-1",logfc.threshold=log2(1.1),assay = "Spatial",slot = 'data')
marker_D <-marker_D[which(marker_D[,1]<=0.05),]
up_D <-rownames(marker_D[which(marker_D[,2]>=0),])
down_D <-rownames(marker_D[which(marker_D[,2]<=0),])
length(up_D)
length(down_D)
up <-intersect(up_A,intersect(up_B,intersect(up_C,up_D)))
down <-intersect(down_A,intersect(down_B,intersect(down_C,down_D)))
#
marker_ABCD <-FindMarkers(brain,ident.1=c("2_A1-2","2_B1-2","2_C1-2","2_D1-2"),ident.2=c("2_A1-1","2_B1-1","2_C1-1","2_D1-1"),logfc.threshold=log2(1.1),assay = "Spatial",slot = 'data')
marker_ABCD <-marker_ABCD[which(marker_ABCD[,1]<=0.05),]
down_ABCD <-rownames(marker_ABCD[which(marker_ABCD[,2]<=0),])
up_ABCD <-rownames(marker_ABCD[which(marker_ABCD[,2]>=0),])
length(up_ABCD)
length(down_ABCD)
intersect(up,up_ABCD)
intersect(down,down_ABCD)
##Venn
setwd("/home/wucheng/Spatial_mouse/Analysis/Merge_data/Different2/C2/Figure_S")
library(venn)
x <- list(up = up_A, up1 = up_B, up2 = up_C,up3 = up_D)
pdf("up.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
x <- list(down = down_A, down1 = down_B, down2 = down_C,down3 = down_D)
pdf("down.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
saveRDS(up,"Up.rds")
saveRDS(down,"Down.rds")

##expression (heatmap plot)
gene <- c(down,up)  
brain@active.ident <-factor(brain$seurat_clusters)
pbmc_c2<-subset(x = brain,idents="2")
pbmc_c2_data <-pbmc_c2[["Spatial"]]@data[gene,]
colnames(pbmc_c2_data) <-pbmc_c2@meta.data[,1]
ind1 <-which(colnames(pbmc_c2_data )=="A1-1")
ind2 <-which(colnames(pbmc_c2_data )=="B1-1")
ind3 <-which(colnames(pbmc_c2_data )=="C1-1")
ind4 <-which(colnames(pbmc_c2_data )=="D1-1")
ind5 <-which(colnames(pbmc_c2_data )=="A1-2")
ind6 <-which(colnames(pbmc_c2_data )=="B1-2")
ind7 <-which(colnames(pbmc_c2_data )=="C1-2")
ind8 <-which(colnames(pbmc_c2_data )=="D1-2")
AA <-pbmc_c2_data[,ind1]
AA1 <-pbmc_c2_data[,ind2]
AA2 <-pbmc_c2_data[,ind3]
AA3 <-pbmc_c2_data[,ind4]
BB <-pbmc_c2_data[,ind5]
BB1 <-pbmc_c2_data[,ind6]
BB2 <-pbmc_c2_data[,ind7]
BB3 <-pbmc_c2_data[,ind8]
mat <-cbind(rowMeans(as.matrix(AA)),rowMeans(as.matrix(AA1)) ,rowMeans(as.matrix(AA2)) ,rowMeans(as.matrix(AA3)),rowMeans(as.matrix(BB)),rowMeans(as.matrix(BB1)) ,rowMeans(as.matrix(BB2)) ,rowMeans(as.matrix(BB3)))
colnames(mat) <-c("A1","B1","C1","D1","A2","B2","C2","D2")
library(pheatmap)
pdf("pheatmap_expression.pdf",width=5,height=15)
pheatmap(mat,cluster_rows = T,cluster_cols=F,show_rownames = T,main = "Heatmap",scale="row")
dev.off()

##Function (GO and KEGG)
dir.create("Function")
setwd("Function")
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
up <-readRDS("/Merge_data/C2/Up.rds")
down <-readRDS("/Merge_data/C2/Down.rds")
gs <- up
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)								
head(ego.bp)
write.csv(ego.bp, file = paste("UP","Go.csv",sep =""))
pdf(paste("UP","Go.pdf",sep =""),width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism   = 'mmu',pvalueCutoff = 0.05)
write.csv(kk,file ="UP_kegg.csv")
pdf("UP_kegg.pdf")
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
gs <- down
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)												
head(ego.bp)
write.csv(ego.bp, file =paste("DOWN","Go.csv",sep =""))
pdf(paste("DOWN","Go.pdf",sep =""),width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism   = 'mmu',pvalueCutoff = 0.05)
write.csv(kk,file ="Down_kegg.csv")
pdf("Down_kegg.pdf")
print(dotplot(kk, showCategory=10,title="KEGG_biological"))
dev.off()

#
# gene set enrichment analysis (GSEA) 
library(msigdbr) #Provide MSigdb database
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) ##mouse
markers <-FindMarkers(brain,ident.1=c("2_A1-2","2_B1-2","2_C1-2","2_D1-2"),ident.2=c("2_A1-1","2_B1-1","2_C1-1","2_D1-1"), min.pct = 0.25, logfc.threshold = 0,assay = "Spatial",slot = 'data')
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) #logFC rank
ranks<- deframe(cluster.genes)
mdb_c2 <- msigdbr(species = "Mus musculus", category = "C2") #C2 gene set 
fgsea_sets = mdb_c2 [grep("^KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name) #select KEGG
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) #run fgsea
saveRDS(fgseaRes,"C2.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") 
pdf('C2_GSEA-fgsea.pdf')
print(p)
dev.off()
