path_data="C:" ###insert path
setwd(path_data)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
n
set.seed(1234567)  # to get same t-SNE.
library(ggplot2)
library(Seurat)
library(pheatmap)
library(DESeq2)
library(umap)
library(Matrix)
library(dplyr)
library(ComplexHeatmap)
library(patchwork)
library(batchelor)

cols <-  c(
  "Blood/dendtritic" = "#b50202",
  "cd163 macrophage" = "#f5b731",
  "dendtritic" = "#db541a",
  "Dendritic/macrophage" = "#f5dd20",
  "PAX5 B cell" = "#db1a95",
  "CD3 T cell/NK cell" = "#c21ddb",
  "Plasmacytoid DC" = "#7d1fe0",
  "cervical_squamousepithelium" = "#1944d1",
  "cervical epithelium KRT7" = "#1944d1",
  "vaginalcervical_squamousepithelium" = "#19b5d1",
  "PECAM1 endothelial/blood" = "#8bd119",
  "Blood.Dendritic1" = "#b50202",
  "Blood.Macrophage1" = "#d2db1a",
  "Blood.Dendritic2" = "#db541a",
  "Blood.Macrophage2" = "#f5dd20",
  "Blood.Macrophage3" = "#8bd119",
  "Blood.Bcell" = "#db1a95",
  "Blood.Tcell" = "#c21ddb",
  "Blood.Dendritic3" = "#f59338",
  "Epi.Cervical" = "#1944d1",
  "Epi.Cervical.KRT7" = "#1944d1",
  "Epi.Vaginal" = "#19b5d1",
  "Endo" = "#7d1fe0",
  "Epi.Cervical1" = "#1944d1",
  "Epi.Cervical2" = "#19b5d1"
)

###load seurat objects

young_seurat10x<-readRDS('5wswabs_ogannot.rds')
old_seurat10x<-readRDS('10wswabs_ogannot.rds')
mid_combinedseurat10x<-readRDS('7combinedwswabs_annot.rds')
plac_seurat10x<-readRDS('plac_annot.rds')

###merge objects, tissue annotations were merged manually for naming consistency

swab.combined <- merge(young_seurat10x, y = c(mid_combinedseurat10x,old_seurat10x), add.cell.ids = c("young", "mid","old"), project = "swabs")
swab.combined <- ScaleData(swab.combined, verbose = FALSE)
swab.combined <- NormalizeData(swab.combined, verbose = FALSE)
swab.combined <- FindVariableFeatures(swab.combined, selection.method = "vst", nfeatures = 10000)
swab.combined <- RunPCA(swab.combined, npcs = 20, verbose = FALSE)
swab.combined <- RunUMAP(swab.combined, reduction = "pca", dims = 1:20)
DimPlot(swab.combined, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)

swab.combined_reclust<-swab.combined
swab.combined_reclust <- FindNeighbors(swab.combined_reclust, dims = 1:10)
swab.combined_reclust <- FindClusters(swab.combined_reclust, resolution = 0.5)
DimPlot(swab.combined_reclust, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)# + NoLegend()

updateidents<-read.csv('newcombined_identrename.csv',stringsAsFactors = F)
Idents(swab.combined_reclust)<-updateidents[match(Idents(swab.combined_reclust),updateidents$orig.ident),2]
saveRDS(swab.combined_reclust, file = "combinedswab_merge_nonorm.rds")
#Identify cell type markers

swabmarkers<-FindAllMarkers(swab.combined_reclust, logfc.threshold = 0.2,min.pct = 0.05,test.use = "wilcox")
write.csv(as.data.frame(swabmarkers), file="newnonorm_combined_markers.csv")

#renaming to retain old identifiers
young_seurat10x_relab<-young_seurat10x
Idents(young_seurat10x_relab)<-unlist(lapply(Idents(young_seurat10x_relab), function(x) paste(x,"_young", sep="")))
old_seurat10x_relab<-old_seurat10x
Idents(old_seurat10x_relab)<-unlist(lapply(Idents(old_seurat10x_relab), function(x) paste(x,"_old", sep="")))
mid_combinedseurat10x_relab<-mid_combinedseurat10x
Idents(mid_combinedseurat10x_relab)<-unlist(lapply(Idents(mid_combinedseurat10x_relab), function(x) paste(x,"_mid", sep="")))

swab.combined_oldident <- merge(young_seurat10x_relab, y = c(mid_combinedseurat10x_relab,old_seurat10x_relab), add.cell.ids = c("young", "mid","old"), project = "swabs")
swab.combined_oldident <- ScaleData(swab.combined_oldident, verbose = FALSE)
swab.combined_oldident <- NormalizeData(swab.combined_oldident, verbose = FALSE)
swab.combined_oldident <- FindVariableFeatures(swab.combined_oldident, selection.method = "vst", nfeatures = 10000)
swab.combined_oldident <- RunPCA(swab.combined_oldident, npcs = 20, verbose = FALSE)
swab.combined_oldident <- RunUMAP(swab.combined_oldident, reduction = "pca", dims = 1:20)

DimPlot(swab.combined_oldident, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1,cols=c("#b50202","#1944d1","#8bd119"))# + NoLegend()
ggsave("combined_norm_coloured.pdf",width =10, height =7)

### epithelial cell subset analysis
epithelialcells <- subset(swab.combined_reclust, idents=c("Epi.Cervical1","Epi.Cervical2"))
epithelialcells <- ScaleData(epithelialcells, verbose = FALSE)
epithelialcells <- NormalizeData(epithelialcells, verbose = FALSE)
epithelialcells <- RunPCA(epithelialcells, npcs = 20, verbose = FALSE)
epithelialcells <- RunUMAP(epithelialcells, reduction = "pca", dims = 1:20)
epithelialcells <- FindNeighbors(epithelialcells, dims = 1:10)
epithelialcells <- FindClusters(epithelialcells, resolution = 0.1)
DimPlot(epithelialcells, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)# + NoLegend()
ggsave("justepithelialcells_nonormmerged_reclust.pdf",width =10, height = 8)

saveRDS(epithelialcells, file = "epi_swabs.rds")

epimarkers<-FindAllMarkers(epithelialcells, logfc.threshold = 0.2,min.pct = 0.05,test.use = "wilcox")
write.csv(as.data.frame(epimarkers), file="epi_newnonorm_markers.csv")

epithelialcells_oldannot<-epithelialcells
rename_idents<-c(rep('young',length(grep("young",colnames(epithelialcells_oldannot)))),rep('mid',length(grep("mid",colnames(epithelialcells_oldannot)))),rep('old',length(grep("old",colnames(epithelialcells_oldannot)))))
Idents(epithelialcells_oldannot)<-rename_idents
DimPlot(epithelialcells_oldannot, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1,cols=c("#b50202","#1944d1","#8bd119"))

epiagemarkers<-FindAllMarkers(epithelialcells_oldannot, logfc.threshold = 0.2,min.pct = 0.05,test.use = "wilcox")
write.csv(as.data.frame(epiagemarkers), file="epi_newnonorm_age__markers.csv")

### generate trophoblast specific gene plots
trophgenes<-read.csv('trophoblast_genes.csv',stringsAsFactors = F)
FeaturePlot(object=swab.combined_reclust, reduction = "umap",features=trophgenes[,1],pt.size=1)
ggsave("trophgenes_combined.pdf",width =15, height = 10)

###placenta analysis
markertest<-FindAllMarkers(plac_dataset, logfc.threshold = 0.5,min.pct = 0.1,test.use = "wilcox")
TE_ALL<-markertest
celltypes<-c("Cytotrophoblast 1",
             "Cytotrophoblast 2",
             "Cytotrophoblast 3",
             "Syncytiotrophoblast",
             "Macrophage",
             "Blood cell 1",
             "Blood cell 2",
             "Blood cell 3",
             "Villous core stroma 1",
             "Villous core stroma 2",
             "Endothelial cell")
heatmap_invivo_genes<-c()
for (i in 1:length(celltypes)){
  temp_cell<-celltypes[i]
  temp_data<-TE_ALL[which(TE_ALL$cluster==temp_cell),]
  temp_data<-temp_data[temp_data$p_val_adj<0.0001,]
  temp_data<-temp_data[order(-temp_data$avg_log2FC),]
  heatmap_invivo_genes<-c(heatmap_invivo_genes,temp_data$gene[1:20])
}
#write.csv(heatmap_invivo_genes,file='placheatmap_genes.csv')
trophheatmapgenes<-heatmap_invivo_genes
#new_order<-levels(plac_dataset)[match(celltypes,levels(plac_dataset))]
plac_dataset@active.ident <- factor(x = plac_dataset@active.ident, levels = celltypes)
plac_dataset<-ScaleData(plac_dataset, features=trophheatmapgenes)

DoHeatmap(object = plac_dataset,features=trophheatmapgenes,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("plac_heatmap.pdf",width =10, height = 10)


