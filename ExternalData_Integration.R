library(Seurat)
library(tidyverse)
library(ggplot2)
library(viridis)

s1<-readRDS("tissue_harmony_WTA_VDJ_filtered.rds")
s1$orig.ident<-"iNKT_Human_WTA_Multi"

###BUGAUT ET AL INTEGRATION
thymus_NKT <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM7667188/")
thymus_Tconv <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM7667131/")
thymus_MAIT <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM7667106/")

thymus_NKT <- CreateSeuratObject(counts = thymus_NKT) %>% UpdateSeuratObject()
thymus_Tconv <- CreateSeuratObject(counts = thymus_Tconv) %>% UpdateSeuratObject()
thymus_MAIT <- CreateSeuratObject(counts = thymus_MAIT) %>% UpdateSeuratObject()

thymus_NKT[["percent.mt"]] <- PercentageFeatureSet(thymus_NKT, pattern = "^MT.")
thymus_Tconv[["percent.mt"]] <- PercentageFeatureSet(thymus_Tconv, pattern = "^MT.")
thymus_MAIT[["percent.mt"]] <- PercentageFeatureSet(thymus_MAIT, pattern = "^MT.")

VlnPlot(thymus_NKT, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(thymus_Tconv, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(thymus_MAIT, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

thymus_NKT[["RNA"]]$data <- NormalizeData(thymus_NKT[["RNA"]]$counts)
thymus_Tconv[["RNA"]]$data <- NormalizeData(thymus_Tconv[["RNA"]]$counts)
thymus_MAIT[["RNA"]]$data <- NormalizeData(thymus_MAIT[["RNA"]]$counts)

thymus_NKT <- SCTransform(thymus_NKT, vars.to.regress = c("nCount_RNA")) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20,force.recalc=TRUE) %>%
  FindClusters(resolution = 0.6)

thymus_Tconv <- SCTransform(thymus_Tconv, vars.to.regress = c("nCount_RNA")) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20,force.recalc=TRUE) %>%
  FindClusters(resolution = 0.6)

thymus_MAIT <- SCTransform(thymus_MAIT, vars.to.regress = c("nCount_RNA")) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20,force.recalc=TRUE) %>%
  FindClusters(resolution = 0.6)

DimPlot(thymus_MAIT, label = TRUE)&coord_equal()
DimPlot(thymus_Tconv, label = TRUE)&coord_equal()
DimPlot(thymus_NKT, label = TRUE)&coord_equal()

thymus_MAIT@meta.data$sample<-"thymus_MAIT"
thymus_Tconv@meta.data$sample<-"thymus_Tconv"
thymus_NKT@meta.data$sample<-"thymus_NKT"

##########################
###THYMUS MAIT ANNOTATION#
##########################

FeaturePlot(thymus_MAIT,features=c("SELL","TBX21","CCL5","NKG7","IL23R"),order=TRUE)

GOI=c("DNTT","NR4A1","EGR2","HIVEP3","CCR9","SATB1","ZBTB16","TBX21","CCL5","NKG7","XCL1","MKI67","RORC","IL23R","CCR7",
      "ISG15")

DotPlot(thymus_MAIT,features=GOI)

library(scCustomize)

a<-Clustered_DotPlot(seurat_object = thymus_MAIT,
                     colors_use_exp =c("paleturquoise1","magenta4"), x_lab_rotate = TRUE,flip=TRUE,#plot_km_elbow = FALSE,
                     exp_color_middle=-1.5,features = GOI)


Idents(thymus_MAIT)<-thymus_MAIT@meta.data$seurat_clusters
thymus_MAIT<-RenameIdents(thymus_MAIT,
                          "0"="MAIT1/17",
                          "1"="Imm.C",
                          "2"="Interm.A/B",
                          "3"="Imm.A",
                          "4"="Interm.A/B",
                          "5"="Cycling",
                          "6"="Interm.C",
                          "7"="Imm.B",
                          "8"="Interm.C",
                          "9"="MAIT1/17",
                          "10"="Interferon",
                          "11"="MAIT1/17")

thymus_MAIT@meta.data$celltype<-Idents(thymus_MAIT)
GOI=c("DNTT","NR4A1","EGR2","HIVEP3","CCR9","SATB1","ZBTB16","TBX21","CCL5","NKG7","XCL1","MKI67","RORC","IL23R","CCR7","SELL",
      "ISG15")

a<-DimPlot(thymus_MAIT, label = TRUE,group.by="celltype")&coord_equal()
b<-Clustered_DotPlot(seurat_object = thymus_MAIT,group.by=c("celltype"),
                     colors_use_exp =c("paleturquoise1","magenta4"), x_lab_rotate = TRUE,flip=TRUE,#plot_km_elbow = FALSE,
                     exp_color_middle=-1.5,features = GOI)

a
b

thymus_MAIT$orig.ident<-"thymus_MAIT"
thymus_MAIT$relabeled_clusters<-thymus_MAIT$celltype


##########################
####THYMUS NKT ANNOTATION#
##########################

DimPlot(thymus_NKT, label = TRUE)&coord_equal()

GOI=c("DNTT","NR4A1","EGR2","HIVEP3","CCR9","SATB1","ZBTB16","TBX21","CCL5","NKG7","XCL1","MKI67","RORC","IL23R","CCR7","SELL","ISG15","CD4")

DotPlot(thymus_NKT,features=GOI,assay="RNA",cluster.idents=TRUE,cols =c("paleturquoise1","magenta4"))&coord_flip()

Clustered_DotPlot(seurat_object = thymus_NKT,group.by=c("celltype"),
                  colors_use_exp =c("paleturquoise1","magenta4"), x_lab_rotate = TRUE,flip=TRUE,#plot_km_elbow = FALSE,
                  features = GOI)

Idents(thymus_NKT)<-thymus_NKT@meta.data$seurat_clusters
thymus_NKT<-RenameIdents(thymus_NKT,
                         "0"="CD4+NKT",
                         "1"="Interm.A",
                         "2"="Imm.B",
                         "3"="NKT1/17",
                         "4"="Imm.A",
                         "5"="Interm.A",
                         "6"="Interm.B",
                         "7"="Cycling",
                         "8"="CD4+NKT",
                         "9"="Mix")
thymus_NKT@meta.data$celltype<-Idents(thymus_NKT)
a<-DimPlot(thymus_NKT, label = TRUE,group.by=c("celltype"))&coord_equal()
b<-DotPlot(thymus_NKT,group.by=c("celltype"),
           features=GOI,assay="SCT",cluster.idents=TRUE,
           cols =c("paleturquoise1","magenta4")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
a|b



thymus_NKT@meta.data$relabeled_clusters<-thymus_NKT@meta.data$celltype
thymus_NKT@meta.data$orig.ident<-"thymus_NKT"


###Integrate All Datasets Together
DimPlot(s1,reduction="harmony.20",group.by="relabeled_clusters")&coord_equal()

all_thymus_MAIT_NKT<-merge(s1,y=c(thymus_MAIT,thymus_NKT))

# run standard anlaysis workflow
all_thymus_MAIT_NKT <- NormalizeData(all_thymus_MAIT_NKT)
all_thymus_MAIT_NKT <- FindVariableFeatures(all_thymus_MAIT_NKT)
all_thymus_MAIT_NKT <- ScaleData(all_thymus_MAIT_NKT,vars.to.regress = "nCount_RNA")
all_thymus_MAIT_NKT <- RunPCA(all_thymus_MAIT_NKT)

all_thymus_MAIT_NKT <- FindNeighbors(all_thymus_MAIT_NKT, dims = 1:30, reduction = "pca")
all_thymus_MAIT_NKT <- FindClusters(all_thymus_MAIT_NKT, resolution = 2, cluster.name = "unintegrated_clusters")

all_thymus_MAIT_NKT <- RunUMAP(all_thymus_MAIT_NKT, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(all_thymus_MAIT_NKT, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))&coord_equal()

DimPlot(all_thymus_MAIT_NKT, reduction = "umap.unintegrated", group.by = c("orig.ident","relabeled_clusters"))&coord_equal()

##Perform Integration

all_thymus_MAIT_NKT <- IntegrateLayers(object = all_thymus_MAIT_NKT, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                       verbose = TRUE)

# re-join layers after integration
all_thymus_MAIT_NKT[["RNA"]] <- JoinLayers(all_thymus_MAIT_NKT[["RNA"]])

all_thymus_MAIT_NKT <- FindNeighbors(all_thymus_MAIT_NKT, reduction = "integrated.cca", dims = 1:30)
all_thymus_MAIT_NKT <- FindClusters(all_thymus_MAIT_NKT, resolution = 1)

all_thymus_MAIT_NKT <- RunUMAP(all_thymus_MAIT_NKT, dims = 1:30, reduction = "integrated.cca")

DimPlot(all_thymus_MAIT_NKT, reduction = "umap", group.by = c("orig.ident", "relabeled_clusters"))

DimPlot(all_thymus_MAIT_NKT, reduction = "umap",split.by="orig.ident", group.by = c("relabeled_clusters"),label=TRUE)&coord_equal()

DimPlot(all_thymus_MAIT_NKT, reduction = "umap",split.by="TRAV_TRAJ_TRBV", group.by = c("orig.ident"),ncol=1,label=TRUE)&coord_equal()



###TRY SCT Integration
# split datasets and process without integration
all_thymus_MAIT_NKT[["RNA"]] <- split(all_thymus_MAIT_NKT[["RNA"]], f = all_thymus_MAIT_NKT$orig.ident)
all_thymus_MAIT_NKT <- SCTransform(all_thymus_MAIT_NKT)
all_thymus_MAIT_NKT <- RunPCA(all_thymus_MAIT_NKT)
all_thymus_MAIT_NKT <- RunUMAP(all_thymus_MAIT_NKT, dims = 1:30)
DimPlot(all_thymus_MAIT_NKT, reduction = "umap", group.by = c("orig.ident", "relabeled_clusters"))

# integrate datasets
all_thymus_MAIT_NKT <- IntegrateLayers(object = all_thymus_MAIT_NKT, method = CCAIntegration, normalization.method = "SCT", verbose = TRUE)
all_thymus_MAIT_NKT <- FindNeighbors(all_thymus_MAIT_NKT, reduction = "integrated.dr", dims = 1:30)
all_thymus_MAIT_NKT <- FindClusters(all_thymus_MAIT_NKT, resolution = 0.6)

all_thymus_MAIT_NKT <- RunUMAP(all_thymus_MAIT_NKT, dims = 1:30, reduction = "integrated.dr")
DimPlot(all_thymus_MAIT_NKT, reduction = "integrated.dr", group.by = c("orig.ident", "relabeled_clusters"))

DimPlot(all_thymus_MAIT_NKT, reduction = "umap",split.by="orig.ident", group.by = c("relabeled_clusters"),label=TRUE)&coord_equal()

saveRDS(all_thymus_MAIT_NKT,"Bugaut_thymus_MAIT_NKT_integration.rds")


############################
####THYMUS TCONV ANNOTATION#
############################

DimPlot(thymus_Tconv,label=TRUE)&coord_equal()
GOI=c("DNTT","NR4A1","EGR2","HIVEP3","CCR9","SATB1","ZBTB16","TBX21","CCL5","NKG7","XCL1","MKI67","RORC","IL23R","CCR7","SELL","ISG15","CD4")

DotPlot(thymus_Tconv,features=GOI,assay="RNA",cluster.idents=TRUE,cols =c("paleturquoise1","magenta4"))&coord_flip()

Clustered_DotPlot(seurat_object = thymus_Tconv,group.by=c("seurat_clusters"),
                  colors_use_exp =c("paleturquoise1","magenta4"), x_lab_rotate = TRUE,flip=TRUE,#plot_km_elbow = FALSE,
                  features = GOI)

Idents(thymus_Tconv)<-thymus_Tconv@meta.data$seurat_clusters
thymus_Tconv<-RenameIdents(thymus_Tconv,
                           "0"="NaiveA",
                           "1"="Imm.B",
                           "2"="NaiveB",
                           "3"="Imm.B",
                           "4"="Imm.C",
                           "5"="Mix_Tconv",
                           "6"="NaiveB",
                           "7"="Imm.A",
                           "8"="NaiveC",
                           "9"="Imm.C")
thymus_Tconv@meta.data$celltype<-Idents(thymus_Tconv)
thymus_Tconv@meta.data$relabeled_clusters<-Idents(thymus_Tconv)

a<-DimPlot(thymus_Tconv, label = TRUE,group.by=c("celltype"))&coord_equal()
b<-DotPlot(thymus_Tconv,group.by=c("celltype"),
           features=GOI,assay="SCT",cluster.idents=TRUE,
           cols =c("paleturquoise1","magenta4")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
a|b

thymus_Tconv@meta.data$orig.ident<-"thymus_Tconv"

###Integrate All Datasets Together - from Buguat Et al.,

all_thymus_MAIT_NKT_Tconv<-merge(s1,y=c(thymus_MAIT,thymus_NKT,thymus_Tconv))

# run standard anlaysis workflow
all_thymus_MAIT_NKT_Tconv <- NormalizeData(all_thymus_MAIT_NKT_Tconv)
all_thymus_MAIT_NKT_Tconv <- FindVariableFeatures(all_thymus_MAIT_NKT_Tconv)
all_thymus_MAIT_NKT_Tconv <- ScaleData(all_thymus_MAIT_NKT_Tconv,vars.to.regress = "nCount_RNA")
all_thymus_MAIT_NKT_Tconv <- RunPCA(all_thymus_MAIT_NKT_Tconv)

all_thymus_MAIT_NKT_Tconv <- FindNeighbors(all_thymus_MAIT_NKT_Tconv, dims = 1:30, reduction = "pca")
all_thymus_MAIT_NKT_Tconv <- FindClusters(all_thymus_MAIT_NKT_Tconv, resolution = 2, cluster.name = "unintegrated_clusters")

all_thymus_MAIT_NKT_Tconv <- RunUMAP(all_thymus_MAIT_NKT_Tconv, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))&coord_equal()

DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap.unintegrated", group.by = c("orig.ident","relabeled_clusters"))&coord_equal()

##Perform Integration

all_thymus_MAIT_NKT_Tconv <- IntegrateLayers(object = all_thymus_MAIT_NKT_Tconv, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                             verbose = TRUE)

# re-join layers after integration
all_thymus_MAIT_NKT_Tconv[["RNA"]] <- JoinLayers(all_thymus_MAIT_NKT_Tconv[["RNA"]])

all_thymus_MAIT_NKT_Tconv <- FindNeighbors(all_thymus_MAIT_NKT_Tconv, reduction = "integrated.cca", dims = 1:30)
all_thymus_MAIT_NKT_Tconv <- FindClusters(all_thymus_MAIT_NKT_Tconv, resolution = 1)

all_thymus_MAIT_NKT_Tconv <- RunUMAP(all_thymus_MAIT_NKT_Tconv, dims = 1:30, reduction = "integrated.cca")

DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap", group.by = c("orig.ident", "relabeled_clusters"))

DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap",split.by="orig.ident", group.by = c("relabeled_clusters"),label=TRUE)&coord_equal()

DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap",split.by="TRAV_TRAJ_TRBV", group.by = c("orig.ident"),ncol=1,label=TRUE)&coord_equal()



### SCT Integration
# split datasets and process without integration
all_thymus_MAIT_NKT_Tconv[["RNA"]] <- split(all_thymus_MAIT_NKT_Tconv[["RNA"]], f = all_thymus_MAIT_NKT_Tconv$orig.ident)
all_thymus_MAIT_NKT_Tconv <- SCTransform(all_thymus_MAIT_NKT_Tconv)
all_thymus_MAIT_NKT_Tconv <- RunPCA(all_thymus_MAIT_NKT_Tconv)
all_thymus_MAIT_NKT_Tconv <- RunUMAP(all_thymus_MAIT_NKT_Tconv, dims = 1:30)
DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap", group.by = c("orig.ident", "relabeled_clusters"))

# integrate datasets
all_thymus_MAIT_NKT_Tconv <- IntegrateLayers(object = all_thymus_MAIT_NKT_Tconv, method = CCAIntegration, normalization.method = "SCT", verbose = TRUE)
all_thymus_MAIT_NKT_Tconv <- FindNeighbors(all_thymus_MAIT_NKT_Tconv, reduction = "integrated.dr", dims = 1:30)
#all_thymus_MAIT_NKT <- FindClusters(all_thymus_MAIT_NKT, resolution = 0.6)

all_thymus_MAIT_NKT_Tconv <- RunUMAP(all_thymus_MAIT_NKT_Tconv, dims = 1:30, reduction = "integrated.dr")
DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "integrated.dr", group.by = c("orig.ident", "relabeled_clusters"))


saveRDS(all_thymus_MAIT_NKT_Tconv,"Bugaut_thymus_MAIT_NKT_Tconv_integration.rds")

all_thymus_MAIT_NKT_Tconv<-readRDS("Bugaut_thymus_MAIT_NKT_Tconv_integration.rds")

table(all_thymus_MAIT_NKT_Tconv@meta.data$sample)
all_thymus_MAIT_NKT_Tconv@meta.data %>% select(sample,percent.mt,nCount_RNA,nFeature_RNA) %>%
  group_by(sample) %>% 
  summarise_at(vars("percent.mt", "nCount_RNA", "nFeature_RNA"), mean) 

table(all_thymus_MAIT_NKT_Tconv@meta.data$sample)

mycols=c('1'="#540d6e",
         '0'="#a1286a",
         '2'="#ee4266",
         '3'="#f78a53",
         '5'="#ffd23f",
         '6'="#62bf41",
         '4'="#2ec4b6",
         '7'="#017365",
         '5'="#ffdc3f",
         '5.1'="#f2ff3f",
         '5.2'="#ffdc3f",
         "Cycling"="#83b692",    
         "MAIT1/17"="#00cecb",   
         "Imm.A"="#ff5e5b",         
         "Imm.C"="#bb8588" ,        
         "Interm.A/B"="#d6ce93" ,   
         "Interm.C"="#ea526f"   ,   
         "Imm.B"="#6930c3"     ,   
         "Interferon"="#001845", 
         "CD4+NKT"="#f9ada0", 
         "NKT1/17"="#5b3758", 
         "Interm.A"="#315100", 
         "Mix"="#ff70a6",     
         "Interm.B"="#d0cfec", 
         "NaiveB"="#83b692",  
         "NaiveA"="#e01e37", 
         "Mix_Tconv"="#7de2d1", 
         "NaiveC"="#dbbca0")  

library(ggrepel)
DimPlot(all_thymus_MAIT_NKT_Tconv, reduction = "umap",split.by="orig.ident", ncol=2,group.by = c("relabeled_clusters"),cols=mycols,label=TRUE,repel=TRUE)&coord_equal()#&geom_text_repel()


####
###########
##Single-cell analysis reveals differences among iNKT cells colonizing peripheral organs and identifies Klf2 as a key gene for iNKT emigration
##Wang et al.,
###https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE130184
##https://pubmed.ncbi.nlm.nih.gov/35915069/
data1 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733520/")
data2 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733521/")
data3 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733522/")
data4 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733523/")
data5 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733524/")
data6 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733525/")
data7 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733526/")
data8 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM3733527/")

data9 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908722/")
data10 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908723/")
data11 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908724/")
data12 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908725/")
data13 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908726/")
data14 <- Read10X(data.dir = "/Users/rjayasin/Desktop/OtherData/GSM4908727/")


Initial_Clustering = function(seurat_object,name){
  seurat_object <- CreateSeuratObject(counts = seurat_object) %>% UpdateSeuratObject()
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT.")
  seurat_object[["RNA"]]$data <- NormalizeData(seurat_object[["RNA"]]$counts)
  seurat_object <- SCTransform(seurat_object, vars.to.regress = c("nCount_RNA")) %>%
    RunPCA(npcs = 20) %>%
    RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(reduction = "pca", dims = 1:20,force.recalc=TRUE) %>%
    FindClusters(resolution = 0.6)
  seurat_object@meta.data$sample<-paste0(name)
  #VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
  #DimPlot(seurat_object, label = TRUE)&coord_equal()
  return(seurat_object)
}

data1=Initial_Clustering(data1,"NKT_mm_Thymus_St0_1")
data2=Initial_Clustering(data2,"NKT_mm_Thymus_St0_2")
data3=Initial_Clustering(data3,"NKT_mm_Thymus_St1_1")
data4=Initial_Clustering(data4,"NKT_mm_Thymus_St1_2")
data5=Initial_Clustering(data5,"NKT_mm_Thymus_St2_1")
data6=Initial_Clustering(data6,"NKT_mm_Thymus_St2_2")
data7=Initial_Clustering(data7,"NKT_mm_Thymus_St3_1")
data8=Initial_Clustering(data8,"NKT_mm_Thymus_St3_2")

data9=Initial_Clustering(data9,"Mouse liver iNKT cell sample 1")
data10=Initial_Clustering(data10,"Mouse liver iNKT cell sample 2")
data11=Initial_Clustering(data11,"Mouse lymph node iNKT cell sample 1")
data12=Initial_Clustering(data12,"Mouse lymph node iNKT cell sample 2")
data13=Initial_Clustering(data13,"Mouse spleen iNKT cell sample 1")
data14=Initial_Clustering(data14,"Mouse spleen iNKT cell sample 2")


DimPlot(data1, label = TRUE)&coord_equal()
DimPlot(data2, label = TRUE)&coord_equal()
DimPlot(data3, label = TRUE)&coord_equal()
DimPlot(data4, label = TRUE)&coord_equal()
DimPlot(data5, label = TRUE)&coord_equal()
DimPlot(data6, label = TRUE)&coord_equal()
DimPlot(data7, label = TRUE)&coord_equal()
DimPlot(data8, label = TRUE)&coord_equal()


Wang_merge<-merge(data1,y=c(data2,data3,data4,data5,data6,data7,data8,
                            data9,data10,data11,data12,data13,data14))
Wang_merge <- SCTransform(Wang_merge, vars.to.regress = c("nCount_RNA")) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20,force.recalc=TRUE) %>%
  FindClusters(resolution = 0.6)

DimPlot(Wang_merge,group.by="seurat_clusters", label = TRUE)&coord_equal()
DimPlot(Wang_merge,group.by="seurat_clusters",split.by=c("sample"),ncol=5, label = TRUE)&coord_equal()

FeaturePlot(Wang_merge,features=c("CENPA","DAPL1","IFIT1","GIMAP7","GZMA","CCL5","JUNB","AQP3","STMN1","HIST1H2AP"),order=TRUE)&coord_equal()

saveRDS(Wang_merge,"/Users/rjayasin/Desktop/OtherData/Wang_merge.rds")



###Integrate All Datasets Together
DimPlot(s1,reduction="harmony.20",group.by="relabeled_clusters")&coord_equal()

#Integrate with our dataset
#overlap genes and merge
common.features <- intersect(rownames(s1), rownames(Wang_merge))
length(common.features)

Wang_merge[['SCT']] <- NULL
s1[['SCT']] <- NULL

DefaultAssay(s1)<-"RNA"
DefaultAssay(Wang_merge)<-"RNA"

all<-merge(s1[common.features, ],y=Wang_merge[common.features, ])

DefaultAssay(all)<-'RNA'
all<-JoinLayers(all)
all[["RNA"]] <- split(all[["RNA"]], f = all$orig.ident)

# pre-process dataset (without integration)
all <- NormalizeData(all)
all <- FindVariableFeatures(all)
all <- ScaleData(all)
all <- RunPCA(all)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all)
all <- RunUMAP(all, dims = 1:20)

DimPlot(all, group.by = c("relabeled_clusters", "sample"))
DimPlot(all, group.by = c("relabeled_clusters", "orig.ident"))

all <- IntegrateLayers(object = all, method = CCAIntegration, 
                       orig.reduction = "pca",
                       new.reduction = "integrated.cca", verbose = FALSE)
all <- FindNeighbors(all, reduction = "integrated.cca", dims = 1:20)
all <- FindClusters(all)

all <- RunUMAP(all, reduction = "integrated.cca", dims = 1:20)
DimPlot(all, group.by = c("orig.ident", "celltype"),reduction="integrated.cca")
DimPlot(all, group.by = c("orig.ident", "relabeled_clusters"),reduction="umap")

DimPlot(all, split.by="orig.ident",group.by = c("relabeled_clusters_broad"),reduction="umap")

saveRDS(all,"Wang_merge_integration.rds")

all<-readRDS("Wang_merge_integration.rds")




