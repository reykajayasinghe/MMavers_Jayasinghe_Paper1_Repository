##Pseudotime analysis

library(tidyverse)

#evalaute trajectory in two sections - first half of trajectory and second half
#Genes that change across all pseudotime

###Correlate Regulon data with pseudotime
##RNA Data
#exp.data <- as.matrix(tissue_harmony_withhla@assays$SCT@scale.data)
##AUC SCT Data - if stored as assay
exp.data <- as.matrix(trajectory1@assays$AUC@scale.data)
#Transpose matrix to combine with metadata
exp.rawdata.tr<-t(exp.data) %>% as.data.frame()
#Extract gene names to do corr analysis
#For RNA data - may want to subset by variable genes
#variable_genes = tissue_harmony_withhla@assays$SCT@var.features %>% as.list()
#For all genes/TFs
GOI = colnames(exp.rawdata.tr) %>% as.list
#Add AUC matrix to psuedotime matrix
exp.rawdata.tr$barcodes<-row.names(exp.rawdata.tr)
row.names(exp.rawdata.tr)<-NULL
#Extract metadata of interest from seurat_object
traj_pseudo<-trajectory1@meta.data %>% select(pseudotime,relabeled_clusters_broad)
traj_pseudo$barcodes<-row.names(traj_pseudo)
row.names(traj_pseudo)<-NULL
#Merge pseudotime with AUC data
res_s1=merge(traj_pseudo,exp.rawdata.tr)
res_s1$pseudotime=as.numeric(as.character(unlist(res_s1$pseudotime)))
#Perform pairwise correlation for all Tfs
all_st=NULL
for (gene in GOI){
  res_s2 = res_s1 %>% select(pseudotime,gene)
  test=cor.test(res_s2$pseudotime,res_s2[[gene]],method='pearson')
  st=cbind(gene,test$estimate,test$p.value)
  all_st=rbind(all_st,st)
}
#Perform FDR correction
all_st=as.data.frame(all_st)
colnames(all_st)=c('Gene','Estimate','P_value')
all_st$FDR=p.adjust(all_st$P_value, method='fdr')
all_st=all_st[order(all_st$FDR),]
#Filter data and save
all_st$Estimate<-as.numeric(as.character(all_st$Estimate))
all_st$FDR<-as.numeric(as.character(all_st$FDR))
sig_only$Estimate<-as.numeric(as.character(sig_only$Estimate))
sig_only$FDR<-as.numeric(as.character(sig_only$FDR))
sig_only = all_st %>% filter(FDR < 0.05) %>% arrange(Estimate)
write.table(sig_only,"alldata_pseudotimecorr_AUC_allgenes_traj1.tsv",sep='\t',quote=F,row.names=F)
nrow(all_st)
nrow(sig_only)
traj1=sig_only

#evalaute trajectory in two sections - first half of trajectory and second half
#Genes that change across all pseudotime

###Correlate Regulon data with pseudotime
##RNA Data
#exp.data <- as.matrix(tissue_harmony_withhla@assays$SCT@scale.data)
##AUC SCT Data - if stored as assay
exp.data <- as.matrix(trajectory2@assays$AUC@scale.data)
#Transpose matrix to combine with metadata
exp.rawdata.tr<-t(exp.data) %>% as.data.frame()
#Extract gene names to do corr analysis
#For RNA data - may want to subset by variable genes
#variable_genes = tissue_harmony_withhla@assays$SCT@var.features %>% as.list()
#For all genes/TFs
GOI = colnames(exp.rawdata.tr) %>% as.list
#Add AUC matrix to psuedotime matrix
exp.rawdata.tr$barcodes<-row.names(exp.rawdata.tr)
row.names(exp.rawdata.tr)<-NULL
#Extract metadata of interest from seurat_object
traj_pseudo<-trajectory2@meta.data %>% select(pseudotime,relabeled_clusters_broad)
traj_pseudo$barcodes<-row.names(traj_pseudo)
row.names(traj_pseudo)<-NULL
#Merge pseudotime with AUC data
res_s1=merge(traj_pseudo,exp.rawdata.tr)
res_s1$pseudotime=as.numeric(as.character(unlist(res_s1$pseudotime)))
#Perform pairwise correlation for all Tfs
all_st=NULL
for (gene in GOI){
  res_s2 = res_s1 %>% select(pseudotime,gene)
  test=cor.test(res_s2$pseudotime,res_s2[[gene]],method='pearson')
  st=cbind(gene,test$estimate,test$p.value)
  all_st=rbind(all_st,st)
}
#Perform FDR correction
all_st=as.data.frame(all_st)
colnames(all_st)=c('Gene','Estimate','P_value')
all_st$FDR=p.adjust(all_st$P_value, method='fdr')
all_st=all_st[order(all_st$FDR),]
#Filter data and save
all_st$Estimate<-as.numeric(as.character(all_st$Estimate))
all_st$FDR<-as.numeric(as.character(all_st$FDR))
sig_only$Estimate<-as.numeric(as.character(sig_only$Estimate))
sig_only$FDR<-as.numeric(as.character(sig_only$FDR))
sig_only = all_st %>% filter(FDR < 0.05) %>% arrange(Estimate)
write.table(sig_only,"alldata_pseudotimecorr_AUC_allgenes_traj2.tsv",sep='\t',quote=F,row.names=F)
nrow(all_st)
nrow(sig_only)
traj2=sig_only

##Main Plots - only plotting top and bottom genes for main figure

plot1=traj1 %>% filter(Estimate > 0.5 | Estimate < -0.3) %>%
  arrange(Estimate) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Gene=factor(Gene, levels=Gene)) %>%   # This trick update the factor levels
  ggplot( aes(x=Gene, y=Estimate)) +
  geom_segment( aes(xend=Gene, yend=0)) +
  geom_point(aes(colour = Estimate < 0),
             show.legend = TRUE) +
  scale_y_continuous(breaks = round(seq(-1, 1, by = 0.2),1))+
  #geom_hline(yintercept = 25, linetype = "dashed") + 
  #geom_hline(yintercept = 50, linetype = "dashed")
  coord_flip() +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  ylab("Correlation with Pseudotime (Trajectory 1)")
plot2=traj2 %>% filter(Estimate > 0.3 | Estimate < -0.1) %>%
  arrange(Estimate) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Gene=factor(Gene, levels=Gene)) %>%   # This trick update the factor levels
  ggplot( aes(x=Gene, y=Estimate)) +
  geom_segment( aes(xend=Gene, yend=0)) +
  geom_point(aes(colour = Estimate < 0),
             show.legend = TRUE) +
  scale_y_continuous(breaks = round(seq(-1, 1, by = 0.2),1))+
  
  #geom_hline(yintercept = 25, linetype = "dashed") + 
  #geom_hline(yintercept = 50, linetype = "dashed")
  coord_flip() +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  ylab("Correlation with Pseudotime (Trajectory 2)")

##For full trajectory remove outlier clusters that exist at end of psuedotime analysis
trajectory3<-subset(tissue_harmony_WTA_filtered,idents=c("4","6"),invert=TRUE)
#evalaute trajectory in two sections - first half of trajectory and second half
#Genes that change across all pseudotime

###Correlate Regulon data with pseudotime
##RNA Data
#exp.data <- as.matrix(tissue_harmony_withhla@assays$SCT@scale.data)
##AUC SCT Data - if stored as assay
exp.data <- as.matrix(trajectory3@assays$AUC@scale.data)
#Transpose matrix to combine with metadata
exp.rawdata.tr<-t(exp.data) %>% as.data.frame()
#Extract gene names to do corr analysis
#For RNA data - may want to subset by variable genes
#variable_genes = tissue_harmony_withhla@assays$SCT@var.features %>% as.list()
#For all genes/TFs
GOI = colnames(exp.rawdata.tr) %>% as.list
#Add AUC matrix to psuedotime matrix
exp.rawdata.tr$barcodes<-row.names(exp.rawdata.tr)
row.names(exp.rawdata.tr)<-NULL
#Extract metadata of interest from seurat_object
traj_pseudo<-trajectory3@meta.data %>% select(pseudotime,relabeled_clusters_subcluster)
traj_pseudo$barcodes<-row.names(traj_pseudo)
row.names(traj_pseudo)<-NULL
#Merge pseudotime with AUC data
res_s1=merge(traj_pseudo,exp.rawdata.tr)
res_s1$pseudotime=as.numeric(as.character(unlist(res_s1$pseudotime)))
#Perform pairwise correlation for all Tfs
all_st=NULL
for (gene in GOI){
  res_s2 = res_s1 %>% select(pseudotime,gene)
  test=cor.test(res_s2$pseudotime,res_s2[[gene]],method='pearson')
  st=cbind(gene,test$estimate,test$p.value)
  all_st=rbind(all_st,st)
}
#Perform FDR correction
all_st=as.data.frame(all_st)
colnames(all_st)=c('Gene','Estimate','P_value')
all_st$FDR=p.adjust(all_st$P_value, method='fdr')
all_st=all_st[order(all_st$FDR),]
#Filter data and save
all_st$Estimate<-as.numeric(as.character(all_st$Estimate))
all_st$FDR<-as.numeric(as.character(all_st$FDR))
sig_only$Estimate<-as.numeric(as.character(sig_only$Estimate))
sig_only$FDR<-as.numeric(as.character(sig_only$FDR))
sig_only = all_st %>% filter(FDR < 0.05) %>% arrange(Estimate)
write.table(sig_only,"alldata_pseudotimecorr_AUC_allgenes_traj3_noclusters4_6_07042024.tsv",sep='\t',quote=F,row.names=F)
nrow(all_st)
nrow(sig_only)
traj3=sig_only

##Plot featureplots of interest
a<-FeaturePlot(tissue_harmony_WTA_filtered,features=c("auc_BCL3"),reduction="harmony.20",min.cutoff="q05", max.cutoff="q95",order=TRUE,pt.size=0.5)&coord_equal()&scale_colour_viridis_c(option="magma")
b<-FeaturePlot(tissue_harmony_WTA_filtered,features=c("auc_STAT1"),reduction="harmony.20",min.cutoff="q05", max.cutoff="q95",order=TRUE,pt.size=0.5)&coord_equal()&scale_colour_viridis_c(option="magma")
c<-FeaturePlot(tissue_harmony_WTA_filtered,features=c("auc_EOMES"),reduction="harmony.20",min.cutoff="q05", max.cutoff="q95",order=TRUE,pt.size=0.5)&coord_equal()&scale_colour_viridis_c(option="magma")
d<-FeaturePlot(tissue_harmony_WTA_filtered,features=c("auc_SOX4"),reduction="harmony.20",min.cutoff="q05", max.cutoff="q95",order=TRUE,pt.size=0.5)&coord_equal()&scale_colour_viridis_c(option="magma")



