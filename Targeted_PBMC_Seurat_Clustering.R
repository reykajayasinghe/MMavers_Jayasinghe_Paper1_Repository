###Note clustering for Targeted PBMC Data completed with Seurat v.4.3.0
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

###PBMC Data Analysis
###########################################################
###########################################################

#read in data
pbmc.data <- read.table("Data/154a_targeted_v1.11/Combined_154a-v1-11-default_DBEC_MolsPerCell.csv", skip=0, sep=",", header=TRUE, row.names=1)
#read in sample tag calls
calls <- read.table("Data/154a_targeted_v1.11/154a-v1-11-default_Sample_Tag_Calls.csv", sep=",", header=TRUE, row.names=1)

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
pbmc <- CreateSeuratObject(counts = t(rnadata), project = "pbmc", meta.data=calls)
#plots of genes and counts
# Visualize QC metrics as a violin plot
#pdf(paste0("QC_",sample,".pdf"),height=8,width=8,useDingbats=FALSE)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#dev.off()
#add in abseq data as separate assay
#create a new assay to store Abseq information
ab_assay <- CreateAssayObject(counts = t(abdata))
#store abseq data in same object with RNA
pbmc[["ADT"]] <- ab_assay
Assays(pbmc)
#set default assay to RNA
DefaultAssay(pbmc) <- "RNA"

Idents(pbmc)<-pbmc@meta.data$Sample_Tag
pbmc2_1<-subset(pbmc,idents=c("Undetermined","Multiplet"),invert=TRUE)

pbmc2 <- SCTransform(pbmc2_1, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
pbmc2 <- RunPCA(pbmc2, npcs = 10)
pbmc2 <- RunUMAP(pbmc2, reduction = "pca", dims = 1:10)
pbmc2 <- FindNeighbors(pbmc2, reduction = "pca", dims = 1:10,force.recalc=TRUE)
pbmc2 <- FindClusters(pbmc2, resolution = 0.5)
pbmc2 <- RunUMAP(pbmc2, dims = 1:10)

##Evaluate Macrophage population
FeaturePlot(pbmc2,features=c("FCER1G", "FCGR3A", "LYN", "IL1B"),order=TRUE)&coord_equal()
DimPlot(pbmc2,group.by="seurat_clusters",label=TRUE)

table(pbmc2@meta.data$seurat_clusters)

#Remove cluster which is a suspect macrophage population
Idents(pbmc2)<-pbmc2@meta.data$seurat_clusters
pbmc2_2<-subset(pbmc2,idents=c("11"),invert=TRUE)

###Alternative normalization procedure
resolution.range <- seq(from = 0, to = 1, by = 0.1)
#https://satijalab.org/seurat/articles/sctransform_vignette.html
pbmc2_3 <- SCTransform(pbmc2_2, vars.to.regress = c("nCount_RNA"),return.only.var.genes = F)
pbmc2_3 <- RunPCA(pbmc2_3, npcs = 10)
pbmc2_3 <- RunUMAP(pbmc2_3, reduction = "pca", dims = 1:10)
pbmc2_3 <- FindNeighbors(pbmc2_3, reduction = "pca", dims = 1:10,force.recalc=TRUE)
pbmc2_3 <- FindClusters(pbmc2_3, resolution = resolution.range)
pbmc2_3 <- RunUMAP(pbmc2_3, dims = 1:10)

DimPlot(pbmc2_3,group.by="Sample_Tag",label=TRUE)

clustree(pbmc2_3)

saveRDS(pbmc2_3,"Targeted_PBMC_object.rds")
