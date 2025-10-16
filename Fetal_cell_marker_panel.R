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

#load datasets
swab.combined_reclust<-readRDS('combinedswab_merge_nonorm.rds')
young_seurat10x<-readRDS('5wswabs_ogannot.rds')
old_seurat10x<-readRDS('10wswabs_ogannot.rds')
mid_combinedseurat10x<-readRDS('7combinedwswabs_annot.rds')
epithelialcells<-readRDS('epi_swabs.rds')
plac_dataset<-readRDS('plac_reannot.rds')

## isolate epithelial cells in swabs and placewnta
troph.plac_dataset<-subset(plac_dataset, idents=c("Cytotrophoblast 1","Cytotrophoblast 2","Cytotrophoblast 3","Syncytiotrophoblast"))

#Use CCa to merge datasets
epi.compare.anchors <- FindIntegrationAnchors(object.list = list(epithelialcells,troph.plac_dataset), dims = 1:20, anchor.features = 10000,k.filter = 100)
epi.compare <- IntegrateData(anchorset = epi.compare.anchors, dims = 1:20)
DefaultAssay(epi.compare) <- "integrated"
epi.compare <- ScaleData(epi.compare, verbose = FALSE)
epi.compare <- FindVariableFeatures(epi.compare, selection.method = "vst", nfeatures = 10000)
epi.compare <- RunPCA(epi.compare, npcs = 20, verbose = FALSE)
epi.compare <- RunUMAP(epi.compare, reduction = "pca", dims = 1:20)
DimPlot(epi.compare, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)
saveRDS(epi.compare,'combined_reannot_epitroph.rds')

cols <-  c(
  "Epi.Cervical1" = "#1944d1",
  "Epi.Cervical2" = "#19b5d1",
  "Cytotrophoblast 1" = "#db1a95",
  "Cytotrophoblast 2" = "#c21ddb",
  "Cytotrophoblast 3" = "#7d1fe0",
  "Syncytiotrophoblast" = "#5133ff"
)
DimPlot(epi.compare, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)
###recluster
epi.compare.reannot<-epi.compare
epi.compare.reannot <- FindNeighbors(epi.compare.reannot, dims = 1:10)
epi.compare.reannot <- FindClusters(epi.compare.reannot, resolution = 0.1)
DimPlot(epi.compare.reannot, reduction = "umap",dims = c(1, 2),label.size = 4,pt.size =1)# + NoLegend()
FeaturePlot(object=epi.compare.reannot, reduction = "umap",features=c("KRT7","GATA3","VGLL1","DAB2"),  min.cutoff = 0,max.cutoff = 3, pt.size=1)
ggsave("troph_swab_porentialmarkers.pdf",width =12, height = 10)

##identify differential markers between trophoblast and cervical epithelium
epimarkers2<-FindMarkers(epi.compare,"Cytotrophoblast 2","Epi.Cervical2", logfc.threshold = 0,test.use = "wilcox")
write.csv(as.data.frame(epimarkers2), file="cervical_troph_compare2_markers.csv")
epimarkers3<-FindMarkers(epi.compare,"Cytotrophoblast 3","Epi.Cervical2", logfc.threshold = 0,test.use = "wilcox")
write.csv(as.data.frame(epimarkers3), file="cervical_troph_compare3_markerss.csv")
reclustpositive<-FindMarkers(epi.compare.reannot,"4", logfc.threshold = 0,test.use = "wilcox")
write.csv(as.data.frame(reclustpositive), file="potentialPOS_markers_onown.csv")
reclustpositivecmp<-FindMarkers(epi.compare.reannot,"4","3", logfc.threshold = 0,test.use = "wilcox")
write.csv(as.data.frame(reclustpositive), file="potentialPOS_markers_cervcompare.csv")

###voclano plot
cs3_cs5markers_wilcoxmut <- epimarkers3
siglevel<-0.000001
fclevel<-1
cs3_cs5markers_wilcoxmut <- mutate(cs3_cs5markers_wilcoxmut, sig=ifelse(cs3_cs5markers_wilcoxmut$avg_log2FC>fclevel& cs3_cs5markers_wilcoxmut$p_val_adj <siglevel | cs3_cs5markers_wilcoxmut$avg_log2FC<(-fclevel)& cs3_cs5markers_wilcoxmut$p_val_adj <siglevel  , yes = "Sig", no = "Nonsig")) 
cs3_cs5markers_wilcoxmut[cs3_cs5markers_wilcoxmut$avg_log2FC>fclevel & cs3_cs5markers_wilcoxmut$p_val_adj <siglevel,6]<- "UP"
cs3_cs5markers_wilcoxmut[cs3_cs5markers_wilcoxmut$avg_log2FC<(-fclevel) & cs3_cs5markers_wilcoxmut$p_val_adj <siglevel,6]<- "down"
cs3_cs5markers_wilcoxmut <- cbind(gene=rownames(cs3_cs5markers_wilcoxmut ), cs3_cs5markers_wilcoxmut ) 

volc = ggplot(cs3_cs5markers_wilcoxmut,aes(avg_log2FC, -log2(p_val)))+
  geom_point(aes(col=sig))+        scale_color_manual(values=c("#19b5d1","grey", "#7d1fe0")) +
  theme(panel.background = element_blank())
volc+geom_text_repel(max.overlaps=15,data=cs3_cs5markers_wilcoxmut[cs3_cs5markers_wilcoxmut$sig=="UP" | cs3_cs5markers_wilcoxmut$sig=="down",], aes(label=gene))

ggsave("trophvscervical.pdf",width =10, height = 10)

### count cells and refine marker list

#reload data and find markers between two largest re annotated clusters in merged placental and swab dataset
combined_epitroph<-readRDS('combined_reannot_epitroph.rds')

#Compare epithelial clusters between placental and swab dataset
fetalmakers<-FindMarkers(combined_epitroph,"4","3",test.use = "wilcox")
#calculate proportion expressed in
fetalmakers$proportion<-fetalmakers$pct.1-fetalmakers$pct.2
#filter for significant markers expressed in a majority of cells
filteredfetalmakers<-fetalmakers[fetalmakers$p_val_adj<0.05,]
filteredfetalmakers<-filteredfetalmakers[filteredfetalmakers$avg_log2FC>1,]
filteredfetalmakers<-filteredfetalmakers[filteredfetalmakers$pct.2<0.015,]
genes<-rownames(filteredfetalmakers)

#Compare epithelial clusters between placental and swab dataset
trophmakers<-FindMarkers(combined_epitroph,"1","3",test.use = "wilcox")
trophmakers$proportion<-trophmakers$pct.1-trophmakers$pct.2
filteredtrophmakers<-trophmakers[trophmakers$p_val_adj<0.001,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$avg_log2FC>1,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$pct.2<0.01,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$proportion>0.25,]
genest3<-rownames(filteredtrophmakers)

#Compare epithelial clusters between placental and swab dataset
trophmakers<-FindMarkers(combined_epitroph,"1","0",test.use = "wilcox")
trophmakers$proportion<-trophmakers$pct.1-trophmakers$pct.2
filteredtrophmakers<-trophmakers[trophmakers$p_val_adj<0.001,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$avg_log2FC>1,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$pct.2<0.01,]
filteredtrophmakers<-filteredtrophmakers[filteredtrophmakers$proportion>0.25,]
genest0<-rownames(filteredtrophmakers)

#merge lists generated above to select for conserved genes which is the fetal cell marker panel
genest_filter<-genest3[genest3 %in% genest0] 
matchedgenes<-genes[genes %in% genest_filter] 
write.csv(matchedgenes, file="fetal_genes_matched.csv")

#heatmap of fetal cell marker panel
DoHeatmap(object = combined_epitroph,features=matchedgenes,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("fetalgenes_epicerv_troph_heatmap.pdf",width =10, height = 10)

#expression of fetal cell marker panel in epi cells
epitrophcells<-subset(combined_epitroph, idents=c("1","0","4"))
DoHeatmap(object = epitrophcells,features=matchedgenes,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))
ggsave("fetalgenes_epicerv_troph_shortened_heatmap.pdf",width =10, height = 10)
epicells<-subset(swab.combined_reclust, idents=c("Epi.Cervical1","Epi.Cervical2"))
DoHeatmap(object = epicells,features=matchedgenes,group.by='ident',group.colors=cols) +
  scale_fill_gradientn(colors = c("#00B0F0", "white", "#FF0000"))

# calculate number of cells expressing markers, using cluster enriched for marker panel
potentialfetal<-subset(combined_epitroph, idents=c("4"))
cells<-colnames(potentialfetal)
youngpercent<-length(grep("young",cells))/ncol(young_seurat10x)*100
midpercent<-length(grep("mid",cells))/ncol(mid_combinedseurat10x)*100
oldpercent<-length(grep("old",cells))/ncol(old_seurat10x)*100

rnaseqcounts<-matrix(nrow=3,ncol=2)
rnaseqcounts[,1]<-c(17,87,50)
rnaseqcounts[,2]<-c(youngpercent,midpercent,oldpercent)
colnames(rnaseqcounts)<-c("counts","perc")
rownames(rnaseqcounts)<-c("young","mid","old")

write.csv(rnaseqcounts, file="fetal_counts_poscell.csv")
