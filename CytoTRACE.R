##CytoTRACE
library(Seurat)
library(Matrix)

#Extract Counts Matrix from Seurat Object
counts_matrix <- tissue_harmony_WTA_filtered[["RNA"]]$counts
write.table(counts_matrix, file=paste0("tissue_harmony_WTA_filtered.txt"), sep="\t", quote = FALSE, row.names = TRUE)

##Run CytoTRACE
library(CytoTRACE)
raw_counts<-read.table(file=paste0("/diskmnt/Datasets/mmy_scratch/Dipersio/mmavers/Analysis_012023/Figures/tissue_harmony_WTA_filtered.txt"),sep="\t")
results <- CytoTRACE(raw_counts, ncores = 20,enableFast = FALSE)
cytoout<-as.data.frame(results$CytoTRACE)
colnames(cytoout)<-c("cytotrace")
##barcodes may need to be modified if they have numbers as barcodes
cytoout$rownames<-rownames(cytoout)
cytoout$rownames<-gsub("X","",as.character(cytoout$rownames))
rownames(cytoout)<-cytoout$rownames
cytoout$rownames<-NULL
#Metadata to add to seurat object
write.table(cytoout, file=paste0("tissue_harmony_WTA_filtered_cytotraceoutput.txt"), sep="\t", quote = FALSE, row.names = TRUE)