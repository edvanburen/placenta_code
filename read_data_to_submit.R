###########################
# This code will read in the data files present on Zenodo
# https://zenodo.org/records/10258020
# and process them to produce the final Seurat object
# Used in analysis
##########################

.libPaths("/n/home06/evb/apps/R_4.1.0_packages_MKL")
library(Matrix)
library(gtools)
library(tidyverse)
library(Seurat)
library(hdf5r)
counts_wd<-"/n/holystore01/LABS/xlin/Lab/evb/data/HaeRyung/deliv_Park_10X_022822_counts/"
samples<-c("Sample_Control_1_F","Sample_Control_2_F","Sample_Control_1_M","Sample_Control_2_M"
           ,"Sample_As_1_F","Sample_As_2_F","Sample_As_1_M","Sample_As_2_M")
j<-0
for(samp in samples){
  j<-j+1
  setwd(paste0(counts_wd,samp,"/outs/"))
  temp<-Read10X_h5("filtered_feature_bc_matrix.h5")
  #temp<-Read10X_h5(paste0(samp,"_filtered_feature_bc_matrix.h5"))
  t1<-unlist(strsplit(samp,split="_"))
  ncells<-ncol(temp)
  meta<-data.frame("id"=rep(samp,ncells),"trt_grp"=rep(t1[2],ncells),"sample_number"=rep(t1[3],ncells),"sex"=rep(t1[4],ncells),stringsAsFactors = FALSE)
  meta$sample_number<-as.numeric(meta$sample_number)
  meta$male<-ifelse(meta$sex=="M",1,0)
  rownames(meta)<-colnames(temp)
  t1<-CreateSeuratObject(counts=temp,assay="RNA",min.cells=10,min.features = 200,meta.data = meta)
  assign(paste0("filtered_feature_bc_matrix_",samp),t1)
  print(j)
}

all_merged<-merge(x=filtered_feature_bc_matrix_Sample_As_1_M,
                  y=c(filtered_feature_bc_matrix_Sample_As_1_F,
                      filtered_feature_bc_matrix_Sample_As_2_F,
                      filtered_feature_bc_matrix_Sample_As_2_M,
                      filtered_feature_bc_matrix_Sample_Control_1_F,
                      filtered_feature_bc_matrix_Sample_Control_1_M,
                      filtered_feature_bc_matrix_Sample_Control_2_F,
                      filtered_feature_bc_matrix_Sample_Control_2_M)
                  ,add.cell.ids<-as.character(rep(1:8)))

# # Compute percent mito ratio

all_merged<-PercentageFeatureSet(all_merged,"^mt-", col.name = "percent_mito")

#filter if total count too high or low
#ncount_RNA, marker genes of weird clusters too
#Cluster that only comes from first sample is enterocyte?
all_merged_filter<-subset(all_merged,subset=percent_mito<2.5 & nFeature_RNA>200 & nFeature_RNA<2500)

all_merged_filter2<-all_merged_filter
all.genes<-rownames(all_merged_filter2)

all_merged_filter2<-NormalizeData(all_merged_filter2,assay="RNA")
all_merged_filter2<-ScaleData(all_merged_filter2,features=all.genes)
all_merged_filter2<-SCTransform(all_merged_filter2,vars.to.regress = "percent_mito")


######################################################

# Integration with Seurat
# Determined to be Necessary after
# evaluating result from above

######################################################


samples<-c("Sample_Control_1_F","Sample_Control_2_F","Sample_Control_1_M","Sample_Control_2_M"
           ,"Sample_As_1_F","Sample_As_2_F","Sample_As_1_M","Sample_As_2_M")
j<-0
counts_wd<-"/n/holystore01/LABS/xlin/Lab/evb/data/HaeRyung/deliv_Park_10X_022822_counts/"
for(samp in samples){
  j<-j+1
  setwd(paste0(counts_wd,samp,"/outs/"))
  temp<-Read10X_h5("filtered_feature_bc_matrix.h5")
  #temp<-Read10X_h5(paste0(samp,"_filtered_feature_bc_matrix.h5"))
  t1<-unlist(strsplit(samp,split="_"))
  ncells<-ncol(temp)
  meta<-data.frame("id"=rep(samp,ncells),"trt_grp"=rep(t1[2],ncells),"sample_number"=rep(t1[3],ncells),"sex"=rep(t1[4],ncells),stringsAsFactors = FALSE)
  meta$sample_number<-as.numeric(meta$sample_number)
  meta$male<-ifelse(meta$sex=="M",1,0)
  rownames(meta)<-colnames(temp)
  t1<-CreateSeuratObject(counts=temp,assay="RNA",min.cells=10,min.features = 200,meta.data = meta)
  t1<-PercentageFeatureSet(t1,pattern="^mt-",col.name = "percent_mito")
  t1<-subset(t1,percent_mito<5 & nFeature_RNA>200 & nFeature_RNA<2500)
  t1<-SCTransform(t1,vars.to.regress = "percent_mito")
  assign(paste0("filtered_feature_bc_matrix_",samp),t1)
  print(j)
}

obj.list<-list(filtered_feature_bc_matrix_Sample_As_1_F,
                 filtered_feature_bc_matrix_Sample_As_1_M,
                     filtered_feature_bc_matrix_Sample_As_2_F,
                     filtered_feature_bc_matrix_Sample_As_2_M,
                     filtered_feature_bc_matrix_Sample_Control_1_F,
                     filtered_feature_bc_matrix_Sample_Control_1_M,
                     filtered_feature_bc_matrix_Sample_Control_2_F,
                     filtered_feature_bc_matrix_Sample_Control_2_M)

features <- SelectIntegrationFeatures(object.list = obj.list,nfeatures = 3000)
obj.list<-PrepSCTIntegration(obj.list,anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj.list
                                  ,normalization.method="SCT"
                                  ,anchor.features = features)

all_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(all_integrated)<-"integrated"
all_integrated<-ScaleData(all_integrated)
all_integrated<-NormalizeData(all_integrated,assay="RNA")

all_integrated<-RunPCA(all_integrated,features=rownames(all_integrated),npcs=100)


  ElbowPlot(all_integrated,ndims = 50)
  DefaultAssay(all_integrated)<-"integrated"
  all_integrated<-FindNeighbors(all_integrated,dims=1:25)
  all_integrated<-FindClusters(all_integrated,resolution = .8)
  all_integrated<-RunUMAP(all_integrated,dims=1:25)
  all_integrated<-RunTSNE(all_integrated,dims=1:25)
  
  all_integrated<-PrepSCTFindMarkers(all_integrated)
  
  assign(paste0("all_markers_integ_0.8_mito_",pm,"_SCT"),FindAllMarkers(all_integrated, min.diff.pct = .1,
                                                                        min.pct = .25,logfc.threshold=0.25,assay="SCT"))
  assign(paste0("all_markers_integ_0.8_mito_",pm,"_RNA"),FindAllMarkers(all_integrated, min.diff.pct = .1,
                                                                        min.pct = .25,logfc.threshold=0.25,assay="RNA"))
  

  # Clusters 0,1,2 contain what is likely
  # maternal uterine tissue
  pdf(paste0("UMAP_plot_integrated_mito_",pm,".pdf")
      ,width=8,height=6)
  DimPlot(all_integrated, reduction = "umap",label=TRUE)
  DimPlot(all_integrated, reduction = "umap",group.by="id")
  DimPlot(all_integrated, reduction = "umap",group.by="trt_grp")

  for(samp in unique(all_integrated@meta.data$id)){
    print(DimPlot(subset(all_integrated,id==samp), reduction = "umap",label=TRUE)+
            ggtitle(samp))
  }
  DimPlot(all_integrated, reduction = "tsne",label=TRUE)
  DimPlot(all_integrated, reduction = "tsne",group.by="id")
  DimPlot(all_integrated, reduction = "tsne",group.by="trt_grp")
  for(samp in unique(all_integrated@meta.data$id)){
    print(DimPlot(subset(all_integrated,id==samp), reduction = "tsne",label=TRUE)+
            ggtitle(samp))
  }
  dev.off()
  out_name<-paste0("analysis_dat_mito_",pm)
  assign(eval(out_name),all_integrated)
  
  #################################
  # Remove clusters 0, 1, and 2
  ######################################

  analysis_dat_mito_5_remove012<-subset(analysis_dat_mito_5,idents = setdiff(Idents(analysis_dat_mito_5),c(0,1,2)))
  
  ################################################
  # Perform clustering
  ################################################
  
  
  analysis_dat_mito_5_remove012<-FindNeighbors(analysis_dat_mito_5_remove012,dims=1:25)
  analysis_dat_mito_5_remove012<-FindClusters(analysis_dat_mito_5_remove012,resolution = .8)
  analysis_dat_mito_5_remove012<-RunUMAP(analysis_dat_mito_5_remove012,dims=1:25)
  analysis_dat_mito_5_remove012<-RunTSNE(analysis_dat_mito_5_remove012,dims=1:25)
  
  
  pm=5
  pdf(paste0("UMAP_plot_integrated_mito_",pm,"_remove_012_100522.pdf")
      ,width=8,height=6)
  DimPlot(analysis_dat_mito_5_remove012, reduction = "umap",label=TRUE)
  DimPlot(analysis_dat_mito_5_remove012, reduction = "umap",group.by="id")
  DimPlot(analysis_dat_mito_5_remove012, reduction = "umap",group.by="trt_grp")
  #DimPlot(all_merged_filter_post_analysis, reduction = "umap",split.by="id")
  for(samp in unique(analysis_dat_mito_5_remove012@meta.data$id)){
    print(DimPlot(subset(analysis_dat_mito_5_remove012,id==samp), reduction = "umap",label=TRUE)+
            ggtitle(samp))
  }
  DimPlot(analysis_dat_mito_5_remove012, reduction = "tsne",label=TRUE)
  DimPlot(analysis_dat_mito_5_remove012, reduction = "tsne",group.by="id")
  DimPlot(analysis_dat_mito_5_remove012, reduction = "tsne",group.by="trt_grp")
  for(samp in unique(analysis_dat_mito_5_remove012@meta.data$id)){
    print(DimPlot(subset(analysis_dat_mito_5_remove012,id==samp), reduction = "tsne",label=TRUE)+
            ggtitle(samp))
  }
  dev.off()
  
  ###########################################################################
  # Cell Type Assignment
  ###########################################################################
  
  final_Seurat_object<-RenameIdents(object=analysis_dat_mito_5_remove012,
                                                   "0" = "Endothelial_1",
                                                   "1" = "Spong_1",
                                                   "2" = "Spong_2",
                                                   "3" = "Tropho_Progenitor_1", 
                                                   "4" = "S-TGC_1", 
                                                   "5" = "B Cell", 
                                                   "6" = "Decidual_1", 
                                                   "7" = "Neutrophil_1", 
                                                   "8" = "GlyT", 
                                                   "9" = "Pericyte", 
                                                   "10" = "Endothelial_2", 
                                                   "11" = "Stromal_1", 
                                                   "12" = "Megakaryocyte", 
                                                   "13" = "Fibroblast", 
                                                   "14" = "Erythroid_Precursor_1",
                                                   "15" = "T_Cell",
                                                   "16" = "Neutrophil_2", 
                                                   "17" = "Macrophage_1", 
                                                   "18" = "Endodermal_Cell", 
                                                   "19" = "Decidual_2", 
                                                   "20" = "Labyr_Tropho_1", 
                                                   "21" = "Yolk_Sack",
                                                   "22" = "Labyr_Tropho_2", 
                                                   "23" = "Erythroid_Precursor_2",
                                                   "24" = "Labyr_Tropho_3",
                                                   "25" = "Erythroid_Precursor_3",
                                                   "26" = "S-TGC_2",
                                                   "27" = "Tropho_Progenitor_2",
                                                   "28" = "Stromal_2",
                                                   "29" = "Pericyte_2",
                                                   "30" = "Dendritic",
                                                   "31" = "Basophil",
                                                   "32" = "Macrophage_2",
                                                   "33" = "NK_Cell_1",
                                                   "34" = "Erythroblast",
                                                   "35" = "NK_Cell_2")
  
  

