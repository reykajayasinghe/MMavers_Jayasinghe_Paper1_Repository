#Creating Loom file for Scenic Run

library(SCopeLoomR)

#Below functions based on:
##https://github.com/ding-lab/PanCan_snATAC_publication/tree/main/Fig3/Run_SCENIC
### FUNCTIONS #####
filter_genes_in_ct <- function(counts, cellInfo) {
  possible.ct <- cellInfo$seurat_clusters %>% unique %>% as.character()
  
  genes.to.keep <- sapply(possible.ct, filter_one_ct, counts = counts)
  
  genes.to.keep <- do.call('c', genes.to.keep)
  genes.to.keep <- unique(genes.to.keep)
  return(genes.to.keep)
}

filter_one_ct <- function( ct, counts){
  cat(paste('filtering genes for', ct, '\n'))
  cells.ct <- rownames(subset(cellInfo, seurat_clusters==ct))
  nCells <- length(cells.ct)
  Cellmin <- 0.01 * nCells
  counts.filtered <- counts[,which(colnames(counts) %in% cells.ct)]
  nonzero <- counts.filtered > 0
  # Sums all TRUE values and returns TRUE if more than Cellmin TRUE values per gene
  keep_genes <- Matrix::rowSums(counts.filtered) >= Cellmin
  keep_genes <- names(keep_genes)[keep_genes]
  return(keep_genes)
}

#this function is copied from SCENIC package
add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

# get the data from the object
cellInfo <- tissue_harmony_WTA_filtered@meta.data %>% dplyr::select(relabeled_clusters)
colnames(cellInfo)<-c("seurat_clusters")
#Seurat V4 command below
#counts <- GetAssayData(object = s1, assay = 'RNA',slot = "counts")
#Seurat V5 command below
counts <- LayerData(tissue_harmony_WTA_filtered, assay = "SCT", layer = "counts")
keep_genes <- filter_genes_in_ct(counts, cellInfo)
keep_genes <- keep_genes[!grepl('^MT-', keep_genes)]

DefaultAssay(tissue_harmony_WTA_filtered) <- 'SCT'
#Seurat v4 command below
#exprMat <- t(FetchData(s1, vars = keep_genes, slot = 'data'))
#Seurat v5 command below
exprMat <- t(FetchData(tissue_harmony_WTA_filtered, vars = keep_genes, layer = "counts"))

# save loom
loom <- build_loom(paste0("tissue_harmony_WTA_filtered_SCT.loom"), dgem=exprMat)
#print(head(get_cell_ids(loom)))
#loom
#loom[["col_attrs"]]
#loom[["row_attrs"]]
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)