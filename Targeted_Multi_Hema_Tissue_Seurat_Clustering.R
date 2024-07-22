###Note clustering for Multi-hematopoetic WTA Targeted Data completed with Seurat v.4.3.0

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

#############WTA TARGETED DATA ANALYSIS 
#read in data
pbmc.data <- read.table("_1_Combined_154b-targeted_RSEC_MolsPerCell.csv", skip=0, sep=",", header=TRUE, row.names=1)

#read in sample tag calls
calls <- read.table("_1_154b-targeted_Sample_Tag_Calls.csv", sep=",", header=TRUE, row.names=1)
#remove multiplets and undetermined
#calls <- calls[calls$Sample_Tag != "Multiplet",]
#remove multiplets from data table
data <- pbmc.data[row.names(pbmc.data) %in% row.names(calls),]
#remove abseqs - WTA data
#rnadata <- data[,c(6:30300)]
#abdata <- data[,c(1:5)]
#emove abseqs - targeted data
rnadata <- data[,c(6:440)]
abdata <- data[,c(1:5)]
#remove hlas from further processing
#hlas <- colnames(rnadata[grep("HLA", colnames(rnadata))])
#rnaminushla <- rnadata[,!(colnames(rnadata) %in% hlas)]
#read in object
pbmc <- CreateSeuratObject(counts = t(rnadata), project = "targeted_wta", meta.data=calls)
#plots of genes and counts
# Visualize QC metrics as a violin plot
#pdf(paste0("QC_",sample,".pdf"),height=8,width=8,useDingbats=FALSE)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"),group.by="Sample_Tag", ncol = 2)
#dev.off()
#add in abseq data as separate assay
#create a new assay to store Abseq information
ab_assay <- CreateAssayObject(counts = t(abdata))
#store abseq data in same object with RNA
pbmc[["ADT"]] <- ab_assay
Assays(pbmc)
#set default assay to RNA
DefaultAssay(pbmc) <- "RNA"

##Add sample tag details

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

#filter out undetermined cells and multiplets
Idents(pbmc)<-pbmc@meta.data$Sample_Tag
pbmc = subset(x = pbmc, idents = c("Undetermined","Multiplet"), invert = TRUE)

pbmc2 <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
pbmc2 <- RunPCA(pbmc2, npcs = 20)
pbmc2 <- RunUMAP(pbmc2, reduction = "pca", dims = 1:20)
pbmc2 <- FindNeighbors(pbmc2, reduction = "pca", dims = 1:20,force.recalc=TRUE)
pbmc2 <- FindClusters(pbmc2, resolution = 0.5)
pbmc2 <- RunUMAP(pbmc2, dims = 1:20)

#Identify Macrophage and remove
FeaturePlot(pbmc2,features=c("FCGR3A","TYROBP","CD74","FCER1G","LYN","IL1B"),order=TRUE)&coord_equal()

DimPlot(pbmc2,group.by="seurat_clusters",label=TRUE)

#Remove cluster which is a suspect macrophage population
Idents(pbmc2)<-pbmc2@meta.data$seurat_clusters
pbmc2_2<-subset(pbmc2,idents=c("9"),invert=TRUE)

###Remove cells that dont match appropriate VDJ for iNKT cells
#Add VDJ Data
vdj_table_assay = read.table("_1_154b-targeted_VDJ_perCell.csv", header = TRUE, row.names = 1, sep = ",",na.strings=c("","NA"))
pbmc2_2 = AddMetaData(object = pbmc2_2, metadata = vdj_table_assay)
# Example of using ifelse() in mutate() to create a new column
pbmc2_2@meta.data$TCR_Alpha_Gamma_V_gene_Dominant_Edited=paste0(pbmc2_2@meta.data$TCR_Alpha_Gamma_V_gene_Dominant)
pbmc2_2@meta.data <- pbmc2_2@meta.data %>%
  mutate(TCR_Alpha_Gamma_V_gene_Dominant_Edited = ifelse(is.na(TCR_Alpha_Gamma_V_gene_Dominant), "TRAV Missing", 
                                                         ifelse(TCR_Alpha_Gamma_V_gene_Dominant == "TRAV10*01", "TRAV10*01", "Other TRAV")))
pbmc2_2@meta.data$TCR_Alpha_Gamma_J_gene_Dominant_Edited=paste0(pbmc2_2@meta.data$TCR_Alpha_Gamma_J_gene_Dominant)
pbmc2_2@meta.data <- pbmc2_2@meta.data %>%
  mutate(TCR_Alpha_Gamma_J_gene_Dominant_Edited = ifelse(is.na(TCR_Alpha_Gamma_J_gene_Dominant), "TRAJ Missing", 
                                                         ifelse(TCR_Alpha_Gamma_J_gene_Dominant == "TRAJ18*01", "TRAJ18*01", "Other TRAJ")))
pbmc2_2@meta.data$TCR_Beta_Delta_V_gene_Dominant_Edited=paste0(pbmc2_2@meta.data$TCR_Beta_Delta_V_gene_Dominant)
pbmc2_2@meta.data <- pbmc2_2@meta.data %>%
  mutate(TCR_Beta_Delta_V_gene_Dominant_Edited = ifelse(is.na(TCR_Beta_Delta_V_gene_Dominant), "TRBV Missing", 
                                                        ifelse(TCR_Beta_Delta_V_gene_Dominant == "TRBV25-1*01", "TRBV25-1*01", "Other TRBV")))
pbmc2_2@meta.data$TRAV_TRAJ_TRBV<-paste0(pbmc2_2@meta.data$TCR_Alpha_Gamma_V_gene_Dominant_Edited,";",pbmc2_2@meta.data$TCR_Alpha_Gamma_J_gene_Dominant_Edited,";",pbmc2_2@meta.data$TCR_Beta_Delta_V_gene_Dominant_Edited)

###REMOVE ANY SINGLE CELL DATA THAT DOES NOT MATCH TRAV/TRAJ/TRBV of INKTs
Idents(pbmc2_2)<-pbmc2_2@meta.data$TRAV_TRAJ_TRBV
pbmc2_3 <- subset(pbmc2_2, idents=c( "Other TRAV;Other TRAJ;Other TRBV","Other TRAV;Other TRAJ;TRBV Missing",
                                     "Other TRAV;Other TRAJ;TRBV25-1*01","Other TRAV;TRAJ18*01;TRBV Missing",
                                     "Other TRAV;TRAJ18*01;TRBV25-1*01","TRAV Missing;TRAJ Missing;Other TRBV",
                                     "TRAV10*01;Other TRAJ;Other TRBV","TRAV10*01;Other TRAJ;TRBV Missing",
                                     "TRAV10*01;Other TRAJ;TRBV25-1*01","TRAV10*01;TRAJ18*01;Other TRBV"),invert=TRUE)


###Alternative normalization procedure
resolution.range <- seq(from = 0, to = 1, by = 0.1)
#https://satijalab.org/seurat/articles/sctransform_vignette.html
pbmc2 <- SCTransform(pbmc2_3, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
pbmc2 <- RunPCA(pbmc2, npcs = 20)
pbmc2 <- RunUMAP(pbmc2, reduction = "pca", dims = 1:20)
pbmc2 <- FindNeighbors(pbmc2, reduction = "pca", dims = 1:20,force.recalc=TRUE)
pbmc2 <- FindClusters(pbmc2, resolution = resolution.range)
pbmc2 <- RunUMAP(pbmc2, dims = 1:20)

###Harmony batch correction
pbmc2 <- RunHarmony(object = pbmc2, group.by.vars = "Sample_Tag", 
                    reduction = "pca",assay.use = "SCT", reduction.save = "harmony")

# Select a range of resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.1)

pbmc2 <- pbmc2 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, force.recalc = T) %>%
  FindClusters(resolution = resolution.range) %>%
  identity()

DimPlot(pbmc2,group.by=c("seurat_clusters"),label=TRUE)&coord_equal()
DimPlot(pbmc2,split.by=c("tissue_type"),group.by=c("seurat_clusters"),ncol=2,label=TRUE)&coord_equal()
pbmc2

##Normalized the ADT Data
DefaultAssay(pbmc2)<-"ADT"
pbmc2 <- NormalizeData(pbmc2, normalization.method = "CLR", margin=1, assay = "ADT") 

saveRDS(pbmc2,"wta_targeted.rds")
