
library(Seurat)
library(Matrix)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(utils)
library(tidyverse)


###ADAPTED FROM
#https://github.com/ding-lab/PanCan_snATAC_publication/tree/main/Fig3/Post_SCENIC

regulons=vector(length=50,mode = "list")
for (iter in 1:50){
wdir='Downstream_Analysis/scenic/data/'
loom=open_loom(paste(wdir,'tissue_harmony_WTA_filtered_SCT_',iter,'_auc_mtx.loom',sep=''),mode='r')

# Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")

regulons[[iter]] <- regulonsToGeneLists(regulons_incidMat)
print(iter)
}

all_tfs=NULL
for (iter in 1:50){
    all_tfs=c(all_tfs,names(regulons[[iter]]))
}
stat=table(all_tfs)
stat=as.data.frame(stat)



stat %>% arrange(-Freq) %>% head()


write.table(stat, paste0('TF_frequency.tsv'),sep='\t',quote=F,row.names=F)

stat<-read.table("TF_frequency.tsv",header=TRUE)

sel_tfs=stat$all_tfs[stat$Freq>=30]

all_st=NULL
for (tf in sel_tfs){
    all_gn=NULL
    for (iter in 1:50){
        if(tf %in% names(regulons[[iter]])){
                g_n=regulons[[iter]][names(regulons[[iter]])==tf][[1]]
                all_gn=c(all_gn,g_n)
                }
        }
        st=as.data.frame(table(all_gn))
        st$TF=tf
        all_st=rbind(all_st,st)
        print(tf)
}
all_st=as.data.frame(all_st)
colnames(all_st)[1]=c('Gene')

write.table(all_st, paste0('Target_frequency_byTF.tsv'),sep='\t',row.names=F,quote=F)

iter=25
loom=open_loom(paste('tissue_harmony_WTA_filtered_SCT_',iter,'_auc_mtx.loom',sep=''),mode='r')


###Make a new regulons-array using updated regulons:
tfs=read.table(paste0('TF_frequency.tsv'),sep='\t',header=T)
tfs_s=tfs[tfs$Freq>=30,]
targ=read.table(paste0('Target_frequency_byTF.tsv'),sep='\t',header=T)
targ_s=targ[targ$Freq>=15,]###May need to adapt


regulons_new=vector(mode='list',length=nrow(tfs_s))
for (i in 1: length(tfs_s$all_tfs)){
    tf=tfs_s$all_tfs[i]
    genes=targ_s$Gene[targ_s$TF==tf]
    regulons_new[[i]]=genes
    names(regulons_new)[[i]]=tf
}


exprMatrix=get_dgem(loom)

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=T, nCores=40)

cells_AUC=AUCell_calcAUC(
  regulons_new,
  cells_rankings,
  nCores = 40,
  normAUC = TRUE,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),
  verbose = TRUE
)
mat=getAUC(cells_AUC)

write.table(mat, paste0("UpdatedByFreq.0.6_TFfreq15regulons_CellsAUC.SCT.tsv"),sep='\t',quote=F,row.names=T)

saveRDS(regulons_new, paste0("UpdatedByFreq.0.6_TFfreq15regulons_CellsAUC.SCT.rds"))
