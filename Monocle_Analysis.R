###Monocle Analysis

library(ggpubr)
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(viridis)

s1<-readRDS("=tissue_harmony_WTA_VDJ_filtered.rds")

Seurat_Object_Diet <- DietSeurat(s1, graphs = "umap")
cds <- as.cell_data_set(Seurat_Object_Diet)

Idents(s1)<-s1@meta.data$relabeled_clusters_broad

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))
# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- s1@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
###MIGHT NEED TO CHANGE BELOW BASED ON YOUR FAVORITE REDUCTION
#cds@int_colData@listData$reducedDims$UMAP <- s1@reductions$umap@cell.embeddings
cds@int_colData@listData$reducedDims$UMAP <- s1@reductions$harmony.20@cell.embeddings


# ...3. Learn trajectory graph ------------------------
#cds <- cluster_cells(cds, resolution=1e-3)
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE,
                   close_loop=TRUE,
                   learn_graph_control=list(minimal_branch_len=10,euclidean_distance_ratio=1.5))

#pdf(paste0(workingdir,"/Learn_trajectory_graph.pdf"), width=4, height=4)
plot_cells(cds,
           color_cells_by = 'relabeled_clusters_broad',
           label_groups_by_cluster = TRUE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

#Visually choose Root Nodes
cds <- order_cells(cds)

a<-plot_cells(cds,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)&coord_equal()


b<-plot_cells(cds,
              color_cells_by = "cytotrace",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=1.5)&coord_equal()&scale_color_viridis_c(option="D")

ggarrange(b, a, ncol = 1, nrow = 2)
saveRDS(cds,"tissue_harmony_WTA_VDJ_filtered_monocle.rds")

pseudotimedata=as.data.frame(pseudotime(cds))
colnames(pseudotimedata) <- c('pseudotime') 
write.table(pseudotimedata,"pseudotime_data.tsv",sep="\t",quote=FALSE,row.names=TRUE)


