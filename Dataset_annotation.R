path_data="C:/Users/dylan/Documents/Wheeler_lab/Fetal_cell_project/RNAseq/"
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

###load datasets from process cellranger output. old: week 10.3, mid: week 7.6, young: week 5.2
old_rawdata<-Read10X(
  'C:/Users/Dylan/Documents/Wheeler lab/Fetal cell project/RNAseq/analysis_4/raw_feature_bc_matrix')
mid_rawdata<-Read10X(
  'C:/Users/Dylan/Documents/Wheeler lab/Fetal cell project/RNAseq/analysis_2/swab/raw_feature_bc_matrix')
young_rawdata<-Read10X(
  'C:/Users/Dylan/Documents/Wheeler lab/Fetal cell project/RNAseq/analysis_3/raw_feature_bc_matrix')
mid2_rawdata<-Read10X(
  'C:/Users/Dylan/Documents/Wheeler lab/Fetal cell project/RNAseq/analysis_1/raw_feature_bc_matrix')
plac_rawdata<-Read10X(
  'C:/Users/Dylan/Documents/Wheeler lab/Fetal cell project/RNAseq/analysis_2/plac/raw_feature_bc_matrix')

### Create seurat obects for each dataset, perform initial normalization, and clustering

young_seurat10x <- CreateSeuratObject(counts = young_rawdata, assay = "RNA", min.cells = 1, min.features = 200)
young_seurat10x[["percent.mt"]] <- PercentageFeatureSet(young_seurat10x, pattern = "^MT-")
VlnPlot(young_seurat10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(young_seurat10x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(young_seurat10x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
young_seurat10x <- subset(young_seurat10x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
young_seurat10x <- NormalizeData(young_seurat10x, verbose = FALSE)
young_seurat10x <- FindVariableFeatures(young_seurat10x, selection.method = "vst", nfeatures = 1000)
young_seurat10x <- ScaleData(young_seurat10x, verbose = FALSE)
young_seurat10x <- RunPCA(young_seurat10x, npcs = 20, verbose = FALSE)
young_seurat10x <- RunUMAP(young_seurat10x, reduction = "pca", dims = 1:20)
young_seurat10x <- FindNeighbors(young_seurat10x, dims = 1:10)
young_seurat10x <- FindClusters(young_seurat10x, resolution = 0.5)

mid_seurat10x <- CreateSeuratObject(counts = mid_rawdata, assay = "RNA", min.cells = 1, min.features = 200)
mid_seurat10x[["percent.mt"]] <- PercentageFeatureSet(mid_seurat10x, pattern = "^MT-")
VlnPlot(mid_seurat10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(mid_seurat10x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mid_seurat10x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
mid_seurat10x <- subset(mid_seurat10x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mid_seurat10x <- NormalizeData(mid_seurat10x, verbose = FALSE)
mid_seurat10x <- FindVariableFeatures(mid_seurat10x, selection.method = "vst", nfeatures = 1000)
mid_seurat10x <- ScaleData(mid_seurat10x, verbose = FALSE)
mid_seurat10x <- RunPCA(mid_seurat10x, npcs = 20, verbose = FALSE)
mid_seurat10x <- RunUMAP(mid_seurat10x, reduction = "pca", dims = 1:20)
mid_seurat10x <- FindNeighbors(mid_seurat10x, dims = 1:10)
mid_seurat10x <- FindClusters(mid_seurat10x, resolution = 0.3)


mid2_seurat10x <- CreateSeuratObject(counts = mid2_rawdata, assay = "RNA", min.cells = 1, min.features = 200)
mid2_seurat10x[["percent.mt"]] <- PercentageFeatureSet(mid2_seurat10x, pattern = "^MT-")
VlnPlot(mid2_seurat10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(mid2_seurat10x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mid2_seurat10x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
mid2_seurat10x <- subset(mid2_seurat10x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mid2_seurat10x <- NormalizeData(mid2_seurat10x, verbose = FALSE)
mid2_seurat10x <- FindVariableFeatures(mid2_seurat10x, selection.method = "vst", nfeatures = 1000)
mid2_seurat10x <- ScaleData(mid2_seurat10x, verbose = FALSE)
mid2_seurat10x <- RunPCA(mid2_seurat10x, npcs = 20, verbose = FALSE)
mid2_seurat10x <- RunUMAP(mid2_seurat10x, reduction = "pca", dims = 1:20)
mid2_seurat10x <- FindNeighbors(mid2_seurat10x, dims = 1:10)
mid2_seurat10x <- FindClusters(mid2_seurat10x, resolution = 0.3)


old_seurat10x <- CreateSeuratObject(counts = old_rawdata, assay = "RNA", min.cells = 1, min.features = 200)
old_seurat10x[["percent.mt"]] <- PercentageFeatureSet(old_seurat10x, pattern = "^MT-")
VlnPlot(old_seurat10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(old_seurat10x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(old_seurat10x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
old_seurat10x <- subset(old_seurat10x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
old_seurat10x <- NormalizeData(old_seurat10x, verbose = FALSE)
old_seurat10x <- FindVariableFeatures(old_seurat10x, selection.method = "vst", nfeatures = 1000)
old_seurat10x <- ScaleData(old_seurat10x, verbose = FALSE)
old_seurat10x <- RunPCA(old_seurat10x, npcs = 20, verbose = FALSE)
old_seurat10x <- RunUMAP(old_seurat10x, reduction = "pca", dims = 1:20)
old_seurat10x <- FindNeighbors(old_seurat10x, dims = 1:10)
old_seurat10x <- FindClusters(old_seurat10x, resolution = 0.5)

plac_seurat10x <- CreateSeuratObject(counts = plac_rawdata, assay = "RNA", min.cells = 1, min.features = 200)
plac_seurat10x[["percent.mt"]] <- PercentageFeatureSet(plac_seurat10x, pattern = "^MT-")
VlnPlot(plac_seurat10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(plac_seurat10x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(plac_seurat10x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plac_seurat10x <- subset(plac_seurat10x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
plac_seurat10x <- NormalizeData(plac_seurat10x, verbose = FALSE)
plac_seurat10x <- FindVariableFeatures(plac_seurat10x, selection.method = "vst", nfeatures = 1000)
plac_seurat10x <- ScaleData(plac_seurat10x, verbose = FALSE)
plac_seurat10x <- RunPCA(plac_seurat10x, npcs = 20, verbose = FALSE)
plac_seurat10x <- RunUMAP(plac_seurat10x, reduction = "pca", dims = 1:20)
plac_seurat10x <- FindNeighbors(plac_seurat10x, dims = 1:10)
plac_seurat10x <- FindClusters(plac_seurat10x, resolution = 0.3)

#combine mid datasets into one

mid.anchors <- FindIntegrationAnchors(object.list = list(mid2_seurat10x,mid_seurat10x), dims = 1:20, anchor.features = 2000,k.filter = 100)
mid.combined <- IntegrateData(anchorset = mid.anchors, dims = 1:20)
DefaultAssay(mid.combined) <- "integrated"
mid.combined <- ScaleData(mid.combined, verbose = FALSE)
mid.combined <- RunPCA(mid.combined, npcs = 20, verbose = FALSE)
mid.combined <- RunUMAP(mid.combined, reduction = "pca", dims = 1:20)
mid.combined <- FindNeighbors(mid.combined, dims = 1:10)
mid.combined <- FindClusters(mid.combined, resolution = 0.3)

#Identify markers for each cluster for annotation. Annotation was performed mannually from these genes
markertest<-FindAllMarkers(young_seurat10x)
write.csv(as.data.frame(markertest), file="young_markers.csv")

markertest<-FindAllMarkers(mid.combined)
write.csv(as.data.frame(markertest), file="midall_markers.csv")

markertest<-FindAllMarkers(old_seurat10x)
write.csv(as.data.frame(markertest), file="old_markers.csv")

markertest<-FindAllMarkers(plac_seurat10x)
write.csv(as.data.frame(markertest), file="plac_markers.csv")

young_seurat10x <- RenameIdents(object = young_seurat10x, `0` = "PECAM1 endothelial/blood")
young_seurat10x <- RenameIdents(object = young_seurat10x, `1` = "PECAM1 endothelial/blood")
young_seurat10x <- RenameIdents(object = young_seurat10x, `2` = "Blood/dendtritic")
young_seurat10x <- RenameIdents(object = young_seurat10x, `3` = "Blood/dendtritic")
young_seurat10x <- RenameIdents(object = young_seurat10x, `4` = "dendtritic")
young_seurat10x <- RenameIdents(object = young_seurat10x, `5` = "cd163 macrophage")
young_seurat10x <- RenameIdents(object = young_seurat10x, `6` = "vaginalcervical_squamousepithelium")
young_seurat10x <- RenameIdents(object = young_seurat10x, `7` = "Dendritic/macrophage")
young_seurat10x <- RenameIdents(object = young_seurat10x, `8` = "CD3 T cell/NK cell")
young_seurat10x <- RenameIdents(object = young_seurat10x, `9` = "PAX5 B cell")
young_seurat10x <- RenameIdents(object = young_seurat10x, `10` = "Plasmacytoid DC")

old_seurat10x <- RenameIdents(object = old_seurat10x, `0` = "PECAM1 endothelial/blood")
old_seurat10x <- RenameIdents(object = old_seurat10x, `1` = "PECAM1 endothelial/blood")
old_seurat10x <- RenameIdents(object = old_seurat10x, `2` = "Blood/dendtritic")
old_seurat10x <- RenameIdents(object = old_seurat10x, `3` = "Dendritic/macrophage")
old_seurat10x <- RenameIdents(object = old_seurat10x, `4` = "Dendritic/macrophage")
old_seurat10x <- RenameIdents(object = old_seurat10x, `5` = "Blood/dendtritic")
old_seurat10x <- RenameIdents(object = old_seurat10x, `6` = "cd163 macrophage")
old_seurat10x <- RenameIdents(object = old_seurat10x, `7` = "vaginalcervical_squamousepithelium")
old_seurat10x <- RenameIdents(object = old_seurat10x, `8` = "CD3 T cell/NK cell")

plac_seurat10x <- RenameIdents(object = plac_seurat10x, `0` = "Red blood cell PTTG1")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `1` = "Cytotrophoblast")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `2` = "Red blood cell")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `3` = "Cytotrophoblast_closetoRBC")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `4` = "PDGFRA stroma mesenchyme")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `5` = "TP63 Troph collumn")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `6` = "cd163 macrophage")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `7` = "CGA high Troph (SYNC?)")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `8` = "PDGFRA neg stroma MMP2")
plac_seurat10x <- RenameIdents(object = plac_seurat10x, `9` = "Endothelial cells")

mid.combined <- RenameIdents(object = mid.combined, `0` = "vaginalcervical_squamousepithelium")
mid.combined <- RenameIdents(object = mid.combined, `1` = "cervical_squamousepithelium")
mid.combined <- RenameIdents(object = mid.combined, `2` = "Blood/dendtritic")
mid.combined <- RenameIdents(object = mid.combined, `3` = "vaginalcervical_squamousepithelium")
mid.combined <- RenameIdents(object = mid.combined, `4` = "Blood/dendtritic")
mid.combined <- RenameIdents(object = mid.combined, `5` = "Blood/dendtritic")
mid.combined <- RenameIdents(object = mid.combined, `6` = "cd163 macrophage")
mid.combined <- RenameIdents(object = mid.combined, `7` = "CD3 T cell/NK cell")
mid.combined <- RenameIdents(object = mid.combined, `8` = "cervical epithelium KRT7")

### save seurat objects for later analysis
saveRDS(plac_seurat10x, file = "plac_annot.rds")
saveRDS(young_seurat10x, file = "5wswabs_ogannot.rds")
saveRDS(old_seurat10x, file = "10wswabs_ogannot.rds")
saveRDS(mid.combined, file = "7combinedwswabs_annot.rds")


