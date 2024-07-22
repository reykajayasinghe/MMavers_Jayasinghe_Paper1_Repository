###Note clustering for Multi-hematopoetic WTA Data completed with Seurat v.4.3.0

set.seed(123)
library(Matrix)
library(Seurat)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(viridis)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)
library(gridExtra)
library(grid)
library(tidyverse)
library(scCustomize)

#read in data
data <- read.table("Combined_WTA-VDJ_RSEC_MolsPerCell.csv", skip=0, sep=",", header=TRUE, row.names=1)
#read in sample tag calls
calls <- read.table("WTA-VDJ_Sample_Tag_Calls.csv", sep=",", header=TRUE, row.names=1)

sample="WTA_version1"

#remove abseqs - WTA data
rnadata <- data[,c(6:32319)]
abdata <- data[,c(1:5)]
#read in object
pbmc <- CreateSeuratObject(counts = t(rnadata), project = "wta", meta.data=calls)
#plots of genes and counts
#Check presence of ZBTB16
rnadata %>% select(ZBTB16)%>% head
#add in abseq data as separate assay
#create a new assay to store Abseq information
ab_assay <- CreateAssayObject(counts = t(abdata))
#store abseq data in same object with RNA
pbmc[["ADT"]] <- ab_assay
Assays(pbmc)
#set default assay to RNA
DefaultAssay(pbmc) <- "RNA"

###Add Meta Data
Idents(pbmc)<-pbmc@meta.data$Sample_Tag
pbmc<-RenameIdents(pbmc,
                   SampleTag01_hs="PB",
                   SampleTag02_hs="PB",
                   SampleTag03_hs="PB",
                   SampleTag04_hs="Cord",
                   SampleTag05_hs="Cord",
                   SampleTag06_hs="Cord",
                   SampleTag07_hs="Thymus",
                   SampleTag08_hs="Thymus",
                   SampleTag09_hs="Thymus",
                   SampleTag10_hs="BM",
                   SampleTag11_hs="BM")
pbmc@meta.data$tissue_type<-Idents(pbmc)

Idents(pbmc)<-pbmc@meta.data$Sample_Tag
pbmc<-RenameIdents(pbmc,
                   SampleTag01_hs="PB1",
                   SampleTag02_hs="PB2",
                   SampleTag03_hs="PB3",
                   SampleTag04_hs="Cord1",
                   SampleTag05_hs="Cord2",
                   SampleTag06_hs="Cord3",
                   SampleTag07_hs="Thymus1",
                   SampleTag08_hs="Thymus2",
                   SampleTag09_hs="Thymus3",
                   SampleTag10_hs="BM1",
                   SampleTag11_hs="BM2")
pbmc@meta.data$tissue_type_donor<-Idents(pbmc)

# store relevant metrics into object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT.", col.name = "percent.mt")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^HB[^(P)]", col.name = "percent.hb")
pbmc <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]", col.name = "percent.RP")

VlnPlot(pbmc,features=c("percent.mt"))+geom_hline(yintercept=20)
#Using 20 % mitochondrial DNA as a filter for cells. 

VlnPlot(pbmc,features=c("nCount_RNA")) + geom_hline(yintercept=20000)

VlnPlot(pbmc,features=c("nFeature_RNA")) + geom_hline(yintercept=4000)

###Filter low quality cells and multiplets
obj = subset(x = pbmc, subset = nFeature_RNA > 200 &
               nFeature_RNA < 4000 &
               nCount_RNA < 20000 & percent.mt<20)

#filter out undetermined cells and multiplets
Idents(obj)<-obj@meta.data$Sample_Tag
obj2 = subset(x = obj, idents = c("Undetermined","Multiplet"), invert = TRUE)

pbmc_withhla <- SCTransform(obj2, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
pbmc_withhla <- RunPCA(pbmc_withhla, npcs = 20)
pbmc_withhla <- RunUMAP(pbmc_withhla, reduction = "pca", dims = 1:20)
pbmc_withhla <- FindNeighbors(pbmc_withhla, reduction = "pca", dims = 1:20,force.recalc=TRUE)
pbmc_withhla <- FindClusters(pbmc_withhla, resolution = 0.5)
pbmc_withhla <- RunUMAP(pbmc_withhla, dims = 1:20)

#Identify Macrophage and remove
FeaturePlot(pbmc_withhla,features=c("FCGR3A","TYROBP","CD74","FCER1G","LYN","IL1B"),order=TRUE)&coord_equal()

DimPlot(pbmc_withhla,group.by="seurat_clusters",label=TRUE)

###Remove Macro/APC Cluster
table(pbmc_withhla@meta.data$seurat_clusters)
Idents(pbmc_withhla)<-pbmc_withhla$seurat_clusters
wta_subset_withhla<-subset(pbmc_withhla,idents=c('13'),invert=TRUE)

###Remove cells that dont match appropriate VDJ for iNKT cells
#Add VDJ Data
vdj_table_assay = read.table("/diskmnt/Datasets/mmy_scratch/Dipersio/mmavers/Analysis_012023/data/WTA-VDJ_VDJ_perCell.csv", header = TRUE, row.names = 1, sep = ",",na.strings=c("","NA"))
wta_subset_withhla = AddMetaData(object = wta_subset_withhla, metadata = vdj_table_assay)
# Example of using ifelse() in mutate() to create a new column
wta_subset_withhla@meta.data$TCR_Alpha_Gamma_V_gene_Dominant_Edited=paste0(wta_subset_withhla@meta.data$TCR_Alpha_Gamma_V_gene_Dominant)
wta_subset_withhla@meta.data <- wta_subset_withhla@meta.data %>%
  mutate(TCR_Alpha_Gamma_V_gene_Dominant_Edited = ifelse(is.na(TCR_Alpha_Gamma_V_gene_Dominant), "TRAV Missing", 
                                                         ifelse(TCR_Alpha_Gamma_V_gene_Dominant == "TRAV10*01", "TRAV10*01", "Other TRAV")))
wta_subset_withhla@meta.data$TCR_Alpha_Gamma_J_gene_Dominant_Edited=paste0(wta_subset_withhla@meta.data$TCR_Alpha_Gamma_J_gene_Dominant)
wta_subset_withhla@meta.data <- wta_subset_withhla@meta.data %>%
  mutate(TCR_Alpha_Gamma_J_gene_Dominant_Edited = ifelse(is.na(TCR_Alpha_Gamma_J_gene_Dominant), "TRAJ Missing", 
                                                         ifelse(TCR_Alpha_Gamma_J_gene_Dominant == "TRAJ18*01", "TRAJ18*01", "Other TRAJ")))
wta_subset_withhla@meta.data$TCR_Beta_Delta_V_gene_Dominant_Edited=paste0(wta_subset_withhla@meta.data$TCR_Beta_Delta_V_gene_Dominant)
wta_subset_withhla@meta.data <- wta_subset_withhla@meta.data %>%
  mutate(TCR_Beta_Delta_V_gene_Dominant_Edited = ifelse(is.na(TCR_Beta_Delta_V_gene_Dominant), "TRBV Missing", 
                                                        ifelse(TCR_Beta_Delta_V_gene_Dominant == "TRBV25-1*01", "TRBV25-1*01", "Other TRBV")))
wta_subset_withhla@meta.data$TRAV_TRAJ_TRBV<-paste0(wta_subset_withhla@meta.data$TCR_Alpha_Gamma_V_gene_Dominant_Edited,";",wta_subset_withhla@meta.data$TCR_Alpha_Gamma_J_gene_Dominant_Edited,";",wta_subset_withhla@meta.data$TCR_Beta_Delta_V_gene_Dominant_Edited)

###REMOVE ANY SINGLE CELL DATA THAT DOES NOT MATCH TRAV/TRAJ/TRBV of INKTs
Idents(wta_subset_withhla)<-wta_subset_withhla@meta.data$TRAV_TRAJ_TRBV
wta_subset_withhla_v2 <- subset(wta_subset_withhla, idents=c( "Other TRAV;Other TRAJ;Other TRBV","Other TRAV;Other TRAJ;TRBV Missing", 
                                                              "Other TRAV;Other TRAJ;TRBV25-1*01","Other TRAV;TRAJ18*01;Other TRBV",
                                                              "Other TRAV;TRAJ18*01;TRBV Missing","Other TRAV;TRAJ18*01;TRBV25-1*01", 
                                                              "TRAV Missing;TRAJ Missing;Other TRBV","TRAV10*01;Other TRAJ;Other TRBV", 
                                                              "TRAV10*01;Other TRAJ;TRBV Missing","TRAV10*01;Other TRAJ;TRBV25-1*01", 
                                                              "TRAV10*01;TRAJ18*01;Other TRBV"),invert=TRUE)

wta_subset_withhla_v2 <- SCTransform(wta_subset_withhla_v2, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
wta_subset_withhla_v2 <- RunPCA(wta_subset_withhla_v2, npcs = 20)
wta_subset_withhla_v2 <- RunUMAP(wta_subset_withhla_v2, reduction = "pca", dims = 1:20)
wta_subset_withhla_v2 <- FindNeighbors(wta_subset_withhla_v2, reduction = "pca", dims = 1:20,force.recalc=TRUE)
wta_subset_withhla_v2 <- FindClusters(wta_subset_withhla_v2, resolution = 0.5)
wta_subset_withhla_v2 <- RunUMAP(wta_subset_withhla_v2, dims = 1:20)

###Harmony batch correction
wta_subset_withhla_v2 <- RunHarmony(object = wta_subset_withhla_v2, group.by.vars = "Sample_Tag", 
                                    reduction = "pca",assay.use = "SCT", reduction.save = "harmony")

# Select a range of resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.1)

wta_subset_withhla_v2 <- wta_subset_withhla_v2 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = resolution.range) %>%
  identity()

DimPlot(wta_subset_withhla_v2,group.by=c("seurat_clusters"),label=TRUE)&coord_equal()
DimPlot(wta_subset_withhla_v2,split.by=c("tissue_type"),group.by=c("seurat_clusters"),ncol=2,label=TRUE)&coord_equal()
wta_subset_withhla_v2

##Normalized the ADT Data
DefaultAssay(wta_subset_withhla_v2)<-"ADT"
wta_subset_withhla_v2 <- NormalizeData(wta_subset_withhla_v2, normalization.method = "CLR", margin=1, assay = "ADT") 

####Run NMF
DefaultAssay(wta_subset_withhla_v2)<-"SCT"
wta_subset_withhla_v2<-wta_subset_withhla_v2 %>% singlet::RunNMF(ncomponents=10)

saveRDS(wta_subset_withhla_v2,"wta_subset_withhla_v2.rds")

