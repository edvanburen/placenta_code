###########################
# This code will produce plots and conduct analysis
# as in the main manuscript
# use read_data_to_submit.R to create the final
# Seurat analysis object or download from
# https://zenodo.org/records/10258020
##########################



# Slight modifications to existing heatmap code
# changes for visualization purposes
plot_heatmap_EVB <- function(dataset, 
                         markers,
                         sort_var = c('seurat_clusters'),
                         n = 8, 
                         anno_var, 
                         anno_colors,
                         hm_limit = c(-2, 0, 2), 
                         hm_colors = c("#4575b4","white","#d73027"),
                         row_font_size = 12) {
  #browser()
  mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
  
  if (is.data.frame(markers)) {
    genes <- get_top_genes(dataset, markers, n)
  } else if (is.character(markers)) {
    genes <- markers
  } else {
    stop('Incorrect input of markers')
  }
  
  
  mat <- mat[match(genes, rownames(mat)),]
  
  anno <- dataset@meta.data %>%
    rownames_to_column(var = "barcode") %>%
    arrange(!!!syms(sort_var))
  
  mat <- t(mat)
  mat <- mat[match(anno$barcode, rownames(mat)),]
  mat <- t(mat)
  
  annos <- list()
  
  for (i in seq_along(1:length(anno_var))) {
    
    err_msg <- paste('Incorrect specification for annotation colors for', anno_var[i])
    
    value <- anno[[anno_var[i]]]
    
    if (is.numeric(value)) {
      
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 'qual'])) {
        
        n <- brewer.pal.info[anno_colors[[i]],]['maxcolors'][[1]]
        pal <- brewer.pal(n = n, name = anno_colors[[i]])
        
        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
                              c(pal[2], pal[(n+1)/2], pal[n-1]))
        
      } else if (length(anno_colors[[i]]) == 3 & all(Scillus:::are_colors(anno_colors[[i]]))) {
        
        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
                              anno_colors[[i]])
      } else {
        stop(err_msg)
      }
      
      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                              col = list(a = col_fun),
                              border = TRUE,
                              annotation_label = anno_var[i],
                              annotation_legend_param=list(labels_gp=gpar(fontsize=row_font_size+4)
                                                           ,title_gp=gpar(fontsize=row_font_size+4
                                                                          ,fontface="bold")
                                                           ,annotation_name_gp=gpar(fontsize=row_font_size+4
                                                                                    ,fontface="bold")))
    } else {
      
      l <- levels(factor(anno[[anno_var[i]]]))
      
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
        
        col <- Scillus:::set_colors(anno_colors[[i]], length(l))
        
      } else if (length(anno_colors[[i]]) >= length(l) & all(Scillus:::are_colors(anno_colors[[i]]))) {
        
        col <- anno_colors[[i]]
        
      } else {
        stop(err_msg)
      }
      
      names(col) <- l
      col <- col[!is.na(names(col))]
      col <- list(a = col)
      
      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                              col = col,
                              border = TRUE,
                              annotation_label = anno_var[i],
                              annotation_name_gp = gpar(fontsize=row_font_size,fontface="bold"),
                              annotation_legend_param=list(labels_gp=gpar(fontsize=row_font_size+4)
                                                           ,title_gp=gpar(fontsize=row_font_size+4
                                                                          ,fontface="bold")
                                                           ,annotation_name_gp=gpar(fontsize=row_font_size+4)))
    }
    names(ha) <- anno_var[i]
    
    annos[[i]] <- ha
  }
  
  annos <- do.call(c, annos)
  
  annos@gap <- rep(unit(1,"mm"), length(annos))
  
  ht <- Heatmap(mat,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(direction = "horizontal",
                                            legend_width = unit(15, "cm"),
                                            title = "Expression",
                                            title_gp=gpar(fontsize=row_font_size),
                                            labels_gp=gpar(fontsize=row_font_size)),
                
                col = colorRamp2(hm_limit, hm_colors),
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = row_font_size),
                top_annotation = annos)
  
  draw(ht, 
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right")
}
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Matrix)
library(gtools)
library(tidyverse)
library(Seurat)
library(metap)
library(ggplot2)
library(MAST)
library(magick)

library(ggpubr)
library(future)
library(Scillus)
load("final_Seurat_object.RData")
dat<-final_Seurat_object


# Before revision the cluster previously
# titled "Hofbauer" was changed to "Macrophage_2"
# to reflect analysis using the Xist gene
# (see below)
# Only run optional code if downloaded object
# from Zenodo
iden<-Idents(final_Seurat_object)
levels(iden)[levels(iden)=="Macrophage"]<-"Macrophage_1"
levels(iden)[levels(iden)=="Hofbauer"]<-"Macrophage_2"
Idents(final_Seurat_object)<-iden
trt_grp<-final_Seurat_object@meta.final_Seurat_objecta$trt_grp

levels(final_Seurat_object@meta.final_Seurat_objecta$cell_type)[levels(final_Seurat_object@meta.final_Seurat_objecta$cell_type)=="Macrophage"]<-"Macrophage_1"

# end optional


dat<-ScaleData(dat,assay="RNA")
dat<-ScaleData(dat,assay="SCT")
DefaultAssay(dat)<-"RNA"

# From a collection of articles and resources
Endothelial_1_marker_genes<-unique(c("Kdr","Pecam1","Plvap","Cdkn1c","Timp3","Cd34","Cldn5"
                              ,"Egfl7","Eng"))
Spong_1_marker_genes<-unique(c("Prl8a9","Flt1","Slco2a1"))
Spong_2_marker_genes<-unique(c("Prl8a9","Flt1"))
Tropho_prog_1_marker_genes<-unique(c("Taf7l","Cited2","Isg20","Rhox6","Rhox9"))
S_TGC_1_marker_genes<-c("Ctsj","Ctsq","Lepr")
B_marker_genes<-c("Igkc","Cd79a")
Decidual_1_marker_genes<-c("Prl8a2","Adm","Cryab")
Neutrophil_1_marker_genes<-c("S100a9","S100a8")
GlyT_marker_genes<-c("Plac8","Pla2g4d","Igfbp7","Ncam1")
Pericyte_1_marker_genes<-c("Acta2","Actg2","Myl9","Mylk")
Endothelial_2_marker_genes<-c("Cdkn1c","Plvap")
Stromal_1_marker_genes<-c("Col4a1")
Megakaryocyte_marker_genes<-c("Pf4","Nrgn","Gp9","Treml1")
Fibroblast_marker_genes<-c("Col1a1","Col3a1","Col6a1","Col5a1","Col5a2","Col6a3")
Erythroid_precursor_1_marker_genes<-c("Hbb-bt","Hba-a1","Hba-a2")
T_marker_genes<-c("Trbc2","Vps37b")
Neutrophil_2_marker_genes<-c("S100a9","S100a8")
Macrophage_1_marker_genes<-c("Apoe","C1qb","Pf4","Ct1a","C1qc")
Endodermal_marker_genes<-c("S100g","Fabp3","Aqp8")
Decidual_2_marker_genes<-c("Prl8a2","Adm","Cryab")
Labyrynthine_tropho_1_marker_genes<-c("Krt8","Krt19")
Yolk_Sack_marker_genes<-c("Apoa4","Fxyd2")
Labyrynthine_tropho_2_marker_genes<-c("Krt8","Krt19","Tfrc")
Erythroid_precursor_2_marker_genes<-c("Hbb-bt","Hba-a1","Hba-a2","Alas2","Bpgm","Snca")
Labyrynthine_tropho_3_marker_genes<-c("Krt8","Krt19")
Erythroid_precursor_3_marker_genes<-c("Hbb-bt","Hba-a1","Hba-a2","Alas2")
S_TGC_2_marker_genes<-c("Ctsj","Ctsq")
Tropho_prog_2_marker_genes<-c("Foxo4","Isg20","Cited2","Taf7l")
Stromal_2_marker_genes<-c("Col4a1")
Pericyte_2_marker_genes<-c("Acta2")
Dendritic_marker_genes<-c("Cd74","Cst3","Cd209a","Cd83")
Basophil_marker_genes<-c("Cd69","Ifitm1","Cd200r3","Mcpt8")
Macrophage_2_marker_genes<-c("Mafb","C1qb","Cd14")
NK_1_marker_genes<-c("Ccl5","Gzma","Nkg7")
Erythroblast_marker_genes<-c("Tfrc","Tmcc2","Slc4a1","Kel","Hmbs")
NK_2_marker_genes<-c("Prf1","Gzmb","Cst7")

all_mark_genes<-unique(do.call(c,lapply(mixedsort(ls(pattern="marker_genes",envir = environment())),get,envir=environment())))

mark_gene_2<-unique(c(Endothelial_1_marker_genes,Spong_1_marker_genes,Spong_2_marker_genes,
               Tropho_prog_1_marker_genes,S_TGC_1_marker_genes,B_marker_genes,
               Decidual_1_marker_genes,Neutrophil_1_marker_genes,GlyT_marker_genes,
               Pericyte_1_marker_genes,Endothelial_2_marker_genes,Stromal_1_marker_genes,
               Megakaryocyte_marker_genes,Fibroblast_marker_genes,Erythroid_precursor_1_marker_genes,
               T_marker_genes,Neutrophil_2_marker_genes,Macrophage_1_marker_genes,
               Endodermal_marker_genes,Decidual_2_marker_genes,Labyrynthine_tropho_1_marker_genes,
               Yolk_Sack_marker_genes,Labyrynthine_tropho_2_marker_genes,
               Erythroid_precursor_2_marker_genes,S_TGC_2_marker_genes,
               Tropho_prog_2_marker_genes,Stromal_2_marker_genes,Pericyte_2_marker_genes,
               Dendritic_marker_genes,Basophil_marker_genes,Macrophage_2_marker_genes,
               NK_1_marker_genes,Erythroblast_marker_genes,NK_2_marker_genes))
nontrop_mark_genes<-unique(c(Endothelial_1_marker_genes,
                            B_marker_genes,
                             Decidual_1_marker_genes,Neutrophil_1_marker_genes,
                             Pericyte_1_marker_genes,Endothelial_2_marker_genes,Stromal_1_marker_genes,
                             Megakaryocyte_marker_genes,Fibroblast_marker_genes,
                            Erythroid_precursor_1_marker_genes,
                             T_marker_genes,Neutrophil_2_marker_genes,Macrophage_1_marker_genes,
                             Endodermal_marker_genes,Decidual_2_marker_genes,
                             Yolk_Sack_marker_genes,
                             Erythroid_precursor_2_marker_genes,Erythroid_precursor_3_marker_genes,
                            Stromal_2_marker_genes,Pericyte_2_marker_genes,
                             Dendritic_marker_genes,Basophil_marker_genes,Macrophage_2_marker_genes,
                             NK_1_marker_genes,Erythroblast_marker_genes,NK_2_marker_genes))

trop_mark_genes<-unique(c(Spong_1_marker_genes,Spong_2_marker_genes,Tropho_prog_1_marker_genes,S_TGC_1_marker_genes,
                          GlyT_marker_genes,Labyrynthine_tropho_1_marker_genes,
                          Labyrynthine_tropho_2_marker_genes,
                          Labyrynthine_tropho_3_marker_genes,S_TGC_2_marker_genes
                          ,Tropho_prog_2_marker_genes))
tropho_idents<-c("Spong_1","Spong_2","Tropho_Progenitor_1","S-TGC_1","GlyT","Labyr_Tropho_1"
                 ,"Labyr_Tropho_2","Labyr_Tropho_3","S-TGC_2","Tropho_Progenitor_2")
nontropho_idents<-setdiff(levels(Idents(dat)),tropho_idents)

############################################################
# Dot Plots
############################################################
dat_subset<-subset(dat,idents=tropho_idents)
dat_inv_subset<-subset(dat,(cell_type%in%c(nontropho_idents)))
pdf("manuscript_dotplots_largefont.pdf"
    ,width=18,height=10)
DotPlot(dat,features=mark_gene_2,assay="RNA")+ theme(axis.text.x=element_text(angle=90, hjust=1),plot.title=element_text(size=rel(1),hjust=.5,face="bold"))+theme(axis.text.x=element_text(angle=90,vjust=.5))+theme(axis.text=element_text(size=20))+theme(axis.title=element_text(size=22,face="bold"),legend.title =element_text(size=18),legend.text = element_text(size=16))+guides(size=guide_legend(title="% Expr."),color=guide_colorbar(title="Avg. Exp."))

DotPlot(dat_subset,features=trop_mark_genes,assay="RNA")+ theme(axis.text.x=element_text(angle=90, hjust=1),plot.title=element_text(size=rel(1),hjust=.5,face="bold"))+theme(axis.text.x=element_text(angle=90,vjust=.5))+theme(axis.text=element_text(size=20))+theme(axis.title=element_text(size=22,face="bold"),legend.title =element_text(size=18),legend.text = element_text(size=16))+guides(size=guide_legend(title="% Expr."),color=guide_colorbar(title="Avg. Exp."))

DotPlot(dat_inv_subset,features=nontrop_mark_genes,assay="RNA")+ theme(axis.text.x=element_text(angle=90, hjust=1),plot.title=element_text(size=rel(1),hjust=.5,face="bold"))+theme(axis.text.x=element_text(angle=90,vjust=.5))+theme(axis.text=element_text(size=16))+theme(axis.title=element_text(size=22,face="bold"),legend.title =element_text(size=18),legend.text = element_text(size=16))+guides(size=guide_legend(title="% Expr."),color=guide_colorbar(title="Avg. Exp."))

dev.off()

############################################################
# Heatmaps
############################################################
dat<-PrepSCTFindMarkers(dat)
#Takes a while
find_all_markers<-FindAllMarkers(dat,logfc.threshold = .1,min.pct=.1,only.pos=TRUE
                                 ,test.use = "MAST",assay="SCT")

tropho_cts<-c("Labyr_Tropho_1","Labyr_Tropho_2","Labyr_Tropho_3","S-TGC_1","S-TGC_2"
              ,"Spong_1","Spong_2","Tropho_Progenitor_1","Tropho_Progenitor_2","GlyT")
nontropho_cts<-setdiff(levels(Idents(dat)),tropho_cts)
DefaultAssay(dat)<-"RNA"
dat_subset<-subset(dat,cell_type%in%c(tropho_cts))
find_all_markers<-find_all_markers
find_all_markers_subset<-find_all_markers%>%filter(cluster%in%c(tropho_cts))
dat_inv_subset<-subset(dat,(cell_type%in%c(nontropho_cts)))
find_all_markers_inv_subset<-find_all_markers%>%filter(cluster%in%c(nontropho_cts))
feat<-c("Prap1","Ctla2a","Guca2b","Gpx3","Apob","Apoa2","Afp")
pdf(file="manuscript_scillus_heatmap_largefont.pdf",width=18,height=14)
plot_heatmap_EVB(dataset = dat,n=3
             ,markers=find_all_markers,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange")))
plot_heatmap_EVB(dataset = dat_subset,n=3
             ,markers=find_all_markers_subset,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange"))
             ,row_font_size = 16)
plot_heatmap_EVB(dataset = dat_inv_subset,n=3
             ,markers=find_all_markers_inv_subset,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange"))
             ,row_font_size = 16)

plot_heatmap_EVB(dataset = dat_subset,n=3
                 ,markers=feat,anno_var = c("cell_type","trt_grp")
                 ,anno_colors=list("Reds",c("blue","orange"))
                 ,row_font_size = 16)
plot_heatmap_EVB(dataset = dat_inv_subset,n=3
                 ,markers=feat,anno_var = c("cell_type","trt_grp")
                 ,anno_colors=list("Reds",c("blue","orange"))
                 ,row_font_size = 16)

dev.off()


# Manually select genes to make sure Prap1 is not added
#seems that it just takes the top 3
genes_manual<-find_all_markers%>%arrange(desc(avg_log2FC))%>%
  group_by(cluster)%>%filter(!gene=="Prap1")%>%filter(row_number()<=3)%>%arrange(cluster)%>%pull(gene)
genes_manual_tropho<-find_all_markers%>%filter(cluster%in%tropho_cts)%>%
  arrange(desc(avg_log2FC))%>%
  group_by(cluster)%>%filter(!gene=="Prap1")%>%filter(row_number()<=3)%>%arrange(cluster)%>%pull(gene)

genes_manual_nontropho<-find_all_markers%>%filter(!cluster%in%tropho_cts)%>%
  arrange(desc(avg_log2FC))%>%
  group_by(cluster)%>%filter(!gene=="Prap1")%>%filter(row_number()<=3)%>%arrange(cluster)%>%pull(gene)

pdf(file="manuscript_heatmap_largefont.pdf",width=18,height=14)
plot_heatmap_EVB(dataset = dat,markers=genes_manual,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange"))
             ,row_font_size = 16)
plot_heatmap_EVB(dataset = dat_subset,markers=genes_manual_tropho,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange"))
             ,row_font_size = 16)
plot_heatmap_EVB(dataset = dat_inv_subset,markers=genes_manual_nontropho,anno_var = c("cell_type","trt_grp")
             ,anno_colors=list("Reds",c("blue","orange"))
             ,row_font_size = 14)
plot_heatmap_EVB(dataset = dat_subset,n=3
                 ,markers=feat,anno_var = c("cell_type","trt_grp")
                 ,anno_colors=list("Reds",c("blue","orange"))
                 ,row_font_size = 16)
plot_heatmap_EVB(dataset = dat_inv_subset,n=3
                 ,markers=feat,anno_var = c("cell_type","trt_grp")
                 ,anno_colors=list("Reds",c("blue","orange"))
                 ,row_font_size = 16)
dev.off()

############################################################
# Barplot
############################################################

#######################
# As vs Control
######################
ct_order<-levels(Idents(dat))
samp<-paste0(dat$sex,"_",dat$trt_grp,"_",dat$sample_number)
ct<-as.character(Idents(dat))
trt_grp<-as.character(dat$trt_grp)
tab_obj<-table(ct)
tab_obj_As<-table(ct[trt_grp=="As"])
tab_obj_Control<-table(ct[trt_grp=="Control"])

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
trt_grp_merge<-data.frame(samp,trt_grp,stringsAsFactors = FALSE)%>%distinct()

df_plot<-data.frame(samp,ct,trt_grp,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot<-left_join(df_plot,trt_grp_merge,by="samp")
df_plot<-df_plot%>%arrange(desc(trt_grp))
df_plot$trt_grp<-factor(df_plot$trt_grp)
df_plot$trt_grp2<-fct_rev(df_plot$trt_grp)

df_plot_orig<-data.frame(samp,ct,trt_grp,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot_orig<-left_join(df_plot_orig,trt_grp_merge,by="samp")
df_plot_orig$trt_grp<-factor(df_plot_orig$trt_grp)
df_plot_orig$trt_grp2<-fct_rev(df_plot_orig$trt_grp)

df_plot_overall<-data.frame(samp,ct,trt_grp,stringsAsFactors = FALSE)
df_plot_overall$trt_grp<-factor(df_plot_overall$trt_grp)
df_plot_overall$trt_grp2<-fct_rev(df_plot_overall$trt_grp)

x_lab<-paste0(ct_order," (N = ",as.numeric(tab_obj[ct_order]),")")
x_lab_by_group<-paste0(ct_order," (N = ",as.numeric(tab_obj_Control[ct_order])
                       ,"; ",as.numeric(tab_obj_As[ct_order]),")")

pdf(file="manuscript_barplot_As_vs_Control.pdf",width=10,height=6)
p1<-ggplot(data=df_plot_overall,aes(x=ct,fill=trt_grp2))+geom_bar(mapping=aes(y=..prop..,group=trt_grp2),position=position_dodge())+geom_point(data=df_plot_orig,mapping=aes(x=ct,y=prop,group=trt_grp2),position=position_dodge(width=0.75))+scale_x_discrete(limits=rev(ct_order),labels=rev(x_lab_by_group),name="Cell Type")+scale_y_continuous(name="Proportion of Cells Per Type Per Mouse (Points are Mouse-Level Averages)",limits = c(0,.15))

df_plot_new<-df_plot_overall%>%group_by(trt_grp)%>%mutate(n_total=n())%>%ungroup()%>%group_by(trt_grp,ct)%>%mutate(n=n())%>%ungroup()%>%distinct(ct,trt_grp,trt_grp2,n_total,n)
p.vals<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_new%>%filter(ct==i)
  p.vals[j]<-prop.test(t1$n,t1$n_total)$p.value
}
names(p.vals)<-ct_order

p.vals2<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_orig%>%filter(ct==i)
  p.vals2[j]<-t.test(t1$prop[t1$trt_grp=="As"],t1$prop[t1$trt_grp=="Control"])$p.value
}
names(p.vals2)<-ct_order

p1+annotate('text',x=1:length(ct_order),y=0.125,label=rev(paste("p =",as.character(formatC(p.vals2,format = "f", digits = 2)))))+coord_flip()+guides(fill=guide_legend(reverse = TRUE))+theme(legend.title = element_blank())

# p1+annotate('text',x=1:length(ct_order),y=0.125,label=rev(paste("p =",as.character(signif(p.vals2,3)))))+coord_flip()+guides(fill=guide_legend(reverse = TRUE))+theme(legend.title = element_blank())
dev.off()

#######################
# M vs. F
######################
ct_order<-levels(Idents(dat))
sex<-dat$sex
samp<-paste0(dat$sex,"_",dat$trt_grp,"_",dat$sample_number)
ct<-as.character(Idents(dat))
trt_grp<-as.character(dat$trt_grp)
tab_obj<-table(ct)
tab_obj_M<-table(ct[sex=="M"])
tab_obj_F<-table(ct[sex=="F"])

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
sex_merge<-data.frame(samp,sex,stringsAsFactors = FALSE)%>%distinct()

df_plot<-data.frame(samp,ct,sex,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot<-left_join(df_plot,sex_merge,by="samp")
df_plot<-df_plot%>%arrange(desc(sex))
df_plot$sex<-factor(df_plot$sex)
df_plot$sex2<-fct_rev(df_plot$sex)

df_plot_orig<-data.frame(samp,ct,sex,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot_orig<-left_join(df_plot_orig,sex_merge,by="samp")
df_plot_orig$sex<-factor(df_plot_orig$sex)
df_plot_orig$sex2<-fct_rev(df_plot_orig$sex)

df_plot_overall<-data.frame(samp,ct,sex,stringsAsFactors = FALSE)
df_plot_overall$sex<-factor(df_plot_overall$sex)
df_plot_overall$sex2<-fct_rev(df_plot_overall$sex)

x_lab<-paste0(ct_order," (N = ",as.numeric(tab_obj[ct_order]),")")
x_lab_by_group<-paste0(ct_order," (N = ",as.numeric(tab_obj_F[ct_order])
                       ,"; ",as.numeric(tab_obj_M[ct_order]),")")

pdf(file="manuscript_barplot_M_vs_F.pdf",width=10,height=6)
p1<-ggplot(data=df_plot_overall,aes(x=ct,fill=sex2))+geom_bar(mapping=aes(y=..prop..,group=sex2),position=position_dodge())+geom_point(data=df_plot_orig,mapping=aes(x=ct,y=prop,group=sex2),position=position_dodge(width=0.75))+scale_x_discrete(limits=rev(ct_order),labels=rev(x_lab_by_group),name="Cell Type")+scale_y_continuous(name="Proportion of Cells Per Type Per Mouse (Points are Mouse-Level Averages)",limits = c(0,.15))

df_plot_new<-df_plot_overall%>%group_by(sex)%>%mutate(n_total=n())%>%ungroup()%>%group_by(sex,ct)%>%mutate(n=n())%>%ungroup()%>%distinct(ct,sex,sex2,n_total,n)
p.vals<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_new%>%filter(ct==i)
  p.vals[j]<-prop.test(t1$n,t1$n_total)$p.value
}
names(p.vals)<-ct_order

p.vals2<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_orig%>%filter(ct==i)
  p.vals2[j]<-t.test(t1$prop[t1$sex=="M"],t1$prop[t1$sex=="F"])$p.value
}
names(p.vals2)<-ct_order
formatC(p.vals2,format = "e", digits = 2)
p1+annotate('text',x=1:length(ct_order),y=0.125,label=rev(paste("p =",as.character(formatC(p.vals2,format = "f", digits = 2)))))+coord_flip()+guides(fill=guide_legend(reverse = TRUE))+theme(legend.title = element_blank())
# p1+annotate('text',x=1:length(ct_order),y=0.125,label=rev(paste("p =",as.character(signif(p.vals2,3)))))+coord_flip()+guides(fill=guide_legend(reverse = TRUE))+theme(legend.title = element_blank())
dev.off()



########################
# Violin Plots
#########################
DefaultAssay(dat)<-"RNA"
impt_feat<-c("Prap1")
ct<-Idents(dat)
dat_AS<-subset(dat,trt_grp=="As")
ct_AS<-Idents(dat_AS)
dat_Control<-subset(dat,trt_grp=="Control")
ct_Control<-Idents(dat_Control)
dat_plot<-subset(dat,features=impt_feat,idents=tropho_idents)
dat_plot2<-subset(dat,features=impt_feat,idents=nontropho_idents)

ct_order<-levels(Idents(dat))

dat_plot3<-subset(dat,features=impt_feat,idents=ct_order[1:17])
dat_plot4<-subset(dat,features=impt_feat,idents=ct_order[18:36])

pdf(file="manuscript_Prap1_expression_As_vs_Control.pdf",width=10,height=6)
for(feat in impt_feat){
  #print(VlnPlot(dat,features=feat,group.by = "sex",split.by="trt_grp",slot="scale.data")+ggtitle(paste0(feat, " Expression in All Cells")))
  #print(VlnPlot(dat_Control,features=feat,group.by = "sex",slot="scale.data")+ggtitle(paste0(feat, " Expression in Control Cells")))
  #print(VlnPlot(dat_AS,features=feat,group.by = "sex",slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells")))
  #print(VlnPlot(dat,features=feat,group.by = "trt_grp",split.by="trt_grp",slot="scale.data")+ggtitle(paste0(feat, " Expression By Treatment and Sex")))
  print(VlnPlot(dat_plot3,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data"))
  print(VlnPlot(dat_plot4,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data"))
  print(VlnPlot(dat_plot,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells From Trophoblast Cell Types")))
  print(VlnPlot(dat_plot2,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells From Non-Trophoblast Cell Types")))
}
dev.off()

pdf(file="manuscript_Prap1_expression_As_cells_M_vs_F.pdf",width=10,height=6)
for(feat in impt_feat){
  #print(VlnPlot(dat,features=feat,group.by = "sex",split.by="trt_grp",slot="scale.data")+ggtitle(paste0(feat, " Expression in All Cells")))
  print(VlnPlot(dat_plot3,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data"))
  print(VlnPlot(dat_plot4,features=feat,group.by = "cell_type",split.by = 'trt_grp',slot="scale.data"))
  print(VlnPlot(dat_Control,features=feat,group.by = "sex",slot="scale.data")+ggtitle(paste0(feat, " Expression in Control Cells")))
  print(VlnPlot(dat_AS,features=feat,group.by = "sex",slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells")))
  print(VlnPlot(dat,features=feat,group.by = "trt_grp",split.by="sex",slot="scale.data")+ggtitle(paste0(feat, " Expression By Treatment and Sex")))
  print(VlnPlot(dat_plot,features=feat,group.by = "cell_type",split.by = 'sex',slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells From Trophoblast Cell Types")))
  print(VlnPlot(dat_plot2,features=feat,group.by = "cell_type",split.by = 'sex',slot="scale.data")+ggtitle(paste0(feat, " Expression in As Cells From Non-Trophoblast Cell Types")))
}
dev.off()

DE_df_out<-all_As_DE_by_ct_MAST%>%filter(p_val_adj<.1)%>%group_by(ct_name)%>%arrange(ct_name,p_val_adj)%>%
  dplyr::rename(pct.express.As.cells=pct.1,pct.express.Control.cells=pct.2)

write_csv(DE_df_out,file="As_vs_Control_DE_100522.csv")



####################################################
#UMAP Plots
#####################################################
dat2<-dat
dat2@meta.data$numeric_ident<-as.numeric(Idents(dat2))
Idents(dat2)<-paste0(dat2@meta.data$numeric_ident,": ",Idents(dat2))
Idents(dat2)<-factor(Idents(dat2),levels=mixedsort(levels(Idents(dat2))))

str_split(paste0(dat2@meta.data$numeric_ident,": ",Idents(dat2)),pattern = ": ")

dat3<-dat
dat3@meta.data$numeric_ident<-as.numeric(Idents(dat3))
Idents(dat3)<-paste0(dat3@meta.data$numeric_ident)
Idents(dat3)<-factor(Idents(dat3),levels=mixedsort(levels(Idents(dat3))))

# try custom
dims = c(1,2);cells = NULL;cols = NULL;pt.size = NULL;reduction = NULL;group.by = NULL;split.by = NULL;shape.by = NULL;order = NULL;shuffle = FALSE;seed = 1;label = FALSE;label.size = 4;label.color = 'black';label.box = FALSE;repel = FALSE;cells.highlight = NULL;cols.highlight = '#DE2D26';sizes.highlight = 1;na.value = 'grey50';ncol = NULL;combine = TRUE;raster = NULL;raster.dpi = c(512, 512)
object=dat;reduction="umap";label=TRUE;object=dat3
DimPlot_EVB <- function(
    object,
    dims = c(1, 2),
    cells = NULL,
    cols = NULL,
    pt.size = NULL,
    reduction = NULL,
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 1,
    label = FALSE,
    label.size = 4,
    label.color = 'black',
    label.box = FALSE,
    repel = FALSE,
    cells.highlight = NULL,
    cols.highlight = '#DE2D26',
    sizes.highlight = 1,
    na.value = 'grey50',
    ncol = NULL,
    combine = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512)
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          labels= NULL,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      plot <- if (is.null(x = orig.groups)) {
        plot + labs(title = NULL)
      } else {
        plot + CenterTitle()
      }
    }
  )
  gsub(".*: ","",levels(Idents(dat2)))
  
  merge_df1<-bind_cols(rownames(dat2@meta.data),dat2@meta.data%>%dplyr::select(numeric_ident,
                                                                     cell_type))
  colnames(merge_df1)[1]<-"cell"
  merge_df1$long_name<-paste0(merge_df1$numeric_ident,": ",merge_df1$cell_type)
  
  merge_df1<-merge_df1%>%dplyr::select(cell,long_name)
  
  plots[[1]]$data$cell<-rownames(plots[[1]]$data)
  plots[[1]]$data<-inner_join(plots[[1]]$data,merge_df1,by="cell")
  plots[[1]]$data$long_name2<-factor(plots[[1]]$data$long_name)
  lev_order<-mixedsort(levels(plots[[1]]$data$long_name2))
  plots[[1]]$data$long_name<-factor(plots[[1]]$data$long_name
                                   ,levels=lev_order)
 
  plots[[1]]<-plots[[1]]+scale_color_discrete(labels=levels(plots[[1]]$data$long_name))
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}
pdf(paste0("manuscript_UMAP_plot.pdf")
    ,width=16,height=8)
DimPlot_EVB(dat3, reduction = "umap",label=TRUE,label.size = 8)
DimPlot(dat, reduction = "umap",group.by="id",label.size = 8)
DimPlot(dat, reduction = "umap",group.by="trt_grp",label.size = 8)
dev.off()

##########################################################################
# Library Complexity
##########################################################################

x_lab<-paste0(ct_order," (N = ",as.numeric(tab_obj[ct_order]),")")
x_lab_by_group<-paste0(ct_order," (N = ",as.numeric(tab_obj_F[ct_order])
                       ,"; ",as.numeric(tab_obj_M[ct_order]),")")

dat2<-dat
nFeature = colSums(x = GetAssayData(object = dat,assay="RNA",slot = "counts") > 0)

pdf(file="library_complexity.pdf",width=14,height=10)
a<-VlnPlot(dat,features="nFeature_RNA",sort=FALSE,assay="RNA",slot="counts"
           ,flip=TRUE,pt.size=0)+NoLegend()+scale_y_continuous(name="Genes Detected Per Cell")+scale_x_discrete(limits=rev(ct_order),labels=rev(x_lab),name="Cell Type")
a<-add_summary(a,fun="median_q1q3")
print(a+coord_flip())

a1<-VlnPlot(dat,features="nFeature_RNA",sort=FALSE
            ,flip=TRUE,pt.size=0,split.by = "trt_grp",split.plot = TRUE)+scale_y_continuous(name="Genes Detected Per Cell",limits=c(0,3200))+scale_x_discrete(limits=rev(ct_order),labels=rev(x_lab_by_group),name="Cell Type")+scale_alpha_continuous(.1)

gr<-a1[1]$data[,"split"]
a1<-add_summary(a1,fun="mean",group = "gr",error.plot="pointrange")
a1<-a1+stat_compare_means(method="t.test",label.y=2900,hide.ns=TRUE,label="p")

print(a1+coord_flip()+guides(fill=guide_legend(reverse = TRUE)))

dev.off()



###########
# Cell Type Proportions
###########

ct_order<-levels(Idents(dat))

tab_obj_As<-table(ct[trt_grp=="As"])
tab_obj_Control<-table(ct[trt_grp=="Control"])

x_lab_by_group<-paste0(ct_order," (N = ",as.numeric(tab_obj_As[ct_order])
                       ,"; ",as.numeric(tab_obj_Control[ct_order]),")")

ct_order_tropho<-ct_order[ct_order%in%tropho_cts]
ct_tropho<-ct[ct%in%tropho_cts]
trt_grp_tropho<-trt_grp[ct%in%tropho_cts]
tab_obj_As_tropho<-table(ct_tropho[trt_grp_tropho=="As"])
tab_obj_Control_tropho<-table(ct_tropho[trt_grp_tropho=="Control"])

x_lab_by_group_tropho<-paste0(ct_order_tropho," (N = ",as.numeric(tab_obj_As_tropho[ct_order_tropho])
                              ,"; ",as.numeric(tab_obj_Control_tropho[ct_order_tropho]),")")

ct_order_nontropho<-ct_order[ct_order%in%nontropho_cts]
ct_nontropho<-ct[ct%in%nontropho_cts]
trt_grp_nontropho<-trt_grp[ct%in%nontropho_cts]
tab_obj_As_nontropho<-table(ct_nontropho[trt_grp_nontropho=="As"])
tab_obj_Control_nontropho<-table(ct_nontropho[trt_grp_nontropho=="Control"])

x_lab_by_group_nontropho<-paste0(ct_order_nontropho," (N = ",as.numeric(tab_obj_As_nontropho[ct_order_nontropho])
                                 ,"; ",as.numeric(tab_obj_Control_nontropho[ct_order_nontropho]),")")

samp<-dat@meta.data$id
trt<-dat@meta.data$trt_grp
ct<-Idents(dat)

gt_merge<-data.frame(samp,trt,stringsAsFactors = FALSE)%>%distinct()

df_plot<-data.frame(samp,ct,trt,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot<-left_join(df_plot,gt_merge,by="samp")
df_plot<-df_plot%>%arrange(desc(trt))
df_plot$trt<-factor(df_plot$trt)
df_plot$trt2<-fct_rev(df_plot$trt)

df_plot_orig<-data.frame(samp,ct,trt,stringsAsFactors = FALSE)%>%group_by(samp,ct)%>%
  summarise(n=n())%>%mutate(prop= n /sum(n))
df_plot_orig<-left_join(df_plot_orig,gt_merge,by="samp")
df_plot_orig$trt<-factor(df_plot_orig$trt)
df_plot_orig$trt2<-fct_rev(df_plot_orig$trt)

df_plot_overall<-data.frame(samp,ct,trt,stringsAsFactors = FALSE)
df_plot_overall$trt<-factor(df_plot_overall$trt)
df_plot_overall$trt2<-fct_rev(df_plot_overall$trt)

ct_order<-levels(Idents(dat))

tab_obj_As<-table(ct[trt_grp=="As"])
tab_obj_Control<-table(ct[trt_grp=="Control"])

x_lab_by_group<-paste0(ct_order," (N = ",as.numeric(tab_obj_As[ct_order])
                       ,"; ",as.numeric(tab_obj_Control[ct_order]),")")

p1<-ggplot(data=df_plot_overall,aes(x=ct,fill=trt2))+geom_bar(mapping=aes(y=..prop..,group=trt2),position=position_dodge())+geom_point(data=df_plot_orig,mapping=aes(x=ct,y=prop,group=trt2),position=position_dodge(width=0.75))+scale_x_discrete(limits=rev(ct_order),labels=rev(x_lab_by_group),name="Cell Type")

df_plot_new<-df_plot_overall%>%group_by(trt)%>%mutate(n_total=n())%>%ungroup()%>%group_by(trt,ct)%>%mutate(n=n())%>%ungroup()%>%distinct(ct,trt,trt2,n_total,n)
p.vals<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_new%>%filter(ct==i)
  p.vals[j]<-prop.test(t1$n,t1$n_total)$p.value
}
names(p.vals)<-ct_order

p.vals2<-rep(NA,length(ct_order))
j<-0
for(i in ct_order){
  j<-j+1
  t1<-df_plot_orig%>%filter(ct==i)
  p.vals2[j]<-t.test(t1$prop[t1$trt=="As"],t1$prop[t1$trt=="Control"])$p.value
}
names(p.vals2)<-ct_order
p.vals2<-format(round(p.vals2,3),nsmall=3)
#p.vals2<-formatC(p.vals2,format='e',digits=2)
library(cowplot)
pdf(file="manuscript_cell_type_proportion.pdf",width=10,height=6)
print(p1+annotate('text',x=1:length(ct_order),y=0.12,label=rev(paste("p =",p.vals2)))+scale_y_continuous(name="Proportion of Cells Per Type Per Mouse (Points are Mouse-Level Averages)")+ylim(0,.13)+theme_cowplot()+coord_flip()+guides(fill=guide_legend(title="Treatment\nGroup",reverse=TRUE)))
dev.off()


############################################################
# Output .csv Files
############################################################
all_As_DE_by_ct_MAST<-tibble()
j<-0
for(ident in levels(Idents(dat))){
  j<-j+1
  print(j)
  c1<-colnames(dat)[trt_grp=="As" & iden==ident]
  c2<-colnames(dat)[trt_grp=="Control" & iden==ident]
  
  dat_subset<-subset(dat,cells = c(c2,c1))
  dat_subset<-PrepSCTFindMarkers(dat_subset)
  t1<-FindMarkers(dat_subset,ident.1="As",assay="SCT"
                  ,group.by="trt_grp",test.use = "MAST",min.pct=0,logfc.threshold = 0)%>%
    mutate(ct_name=ident,comparison="As_vs_Control")
  t1$gene<-rownames(t1)
  all_As_DE_by_ct_MAST<-bind_rows(all_As_DE_by_ct_MAST,t1)
}

all_As_DE_by_ct_MAST$p_val_adj<-p.adjust(all_As_DE_by_ct_MAST$p_val,method="fdr")

# One file
DE_df_out<-all_As_DE_by_ct_MAST%>%filter(p_val_adj<.1)%>%group_by(ct_name)%>%arrange(ct_name,p_val_adj)%>%
  dplyr::rename(pct.express.As.cells=pct.1,pct.express.Control.cells=pct.2)



all_Control_DE_sex_by_ct_MAST<-tibble()
sex<-dat@meta.data$sex
trt_grp<-dat@meta.data$trt_grp
iden<-Idents(dat)
j<-0
for(ident in levels(Idents(dat))){
  j<-j+1
  print(j)
  c1<-colnames(dat)[trt_grp=="Control"]
  
  dat_subset<-subset(dat,cells = c(c1))
  dat_subset<-SCTransform(dat_subset,vars.to.regress = "percent_mito")
  #dat_subset<-PrepSCTFindMarkers(dat_subset)
  #dat_subset<-subset(dat_subset,dat_subset@meta.data$trt_grp=="As")
  t1<-FindMarkers(dat_subset,ident.1="M",assay="SCT"
                  ,group.by="sex",subset.ident=ident,test.use = "MAST",min.pct=0,logfc.threshold = 0)%>%
    mutate(ct_name=ident,comparison="Control_M_vs_F")
  t1$gene<-rownames(t1)
  all_Control_DE_sex_by_ct_MAST<-bind_rows(all_Control_DE_sex_by_ct_MAST,t1)
}

all_Control_DE_sex_by_ct_MAST$p_val_adj<-p.adjust(all_Control_DE_sex_by_ct_MAST$p_val,"fdr")
# another
DE_df_out<-all_Control_DE_sex_by_ct_MAST%>%filter(p_val_adj<.1)%>%group_by(ct_name)%>%arrange(ct_name,p_val_adj)%>%dplyr::rename(pct.express.M.cells=pct.1,pct.express.F.cells=pct.2)

all_As_DE_sex_by_ct_MAST<-tibble()
sex<-dat@meta.data$sex
trt_grp<-dat@meta.data$trt_grp
iden<-Idents(dat)
j<-0
for(ident in levels(Idents(dat))){
  j<-j+1
  print(j)
  c1<-colnames(dat)[trt_grp=="As"]
  
  dat_subset<-subset(dat,cells = c(c1))
  dat_subset<-SCTransform(dat_subset,vars.to.regress = "percent_mito")
  #dat_subset<-PrepSCTFindMarkers(dat_subset)
  #dat_subset<-subset(dat_subset,dat_subset@meta.data$trt_grp=="As")
  t1<-FindMarkers(dat_subset,ident.1="M",assay="SCT"
                  ,group.by="sex",subset.ident=ident,test.use = "MAST",min.pct=0,logfc.threshold = 0)%>%
    mutate(ct_name=ident,comparison="As_M_vs_F")
  t1$gene<-rownames(t1)
  all_As_DE_sex_by_ct_MAST<-bind_rows(all_As_DE_sex_by_ct_MAST,t1)
}

all_As_DE_sex_by_ct_MAST$p_val_adj<-p.adjust(all_As_DE_sex_by_ct_MAST$p_val,"fdr")

DE_df_out<-all_As_DE_sex_by_ct_MAST%>%filter(p_val_adj<.1)%>%group_by(ct_name)%>%arrange(ct_name,p_val_adj)%>%dplyr::rename(pct.express.M.cells=pct.1,pct.express.F.cells=pct.2)
write_csv(DE_df_out,file="manuscript_As_DE_M_vs_F.csv")

# Marker Genes

all_markers_mito_5_remove012_SCT_MAST<-FindAllMarkers(dat, min.diff.pct = .1,test.use="MAST",
                    min.pct = .25,logfc.threshold=0.25,assay="SCT")

manuscript_marker_genes<-all_markers_mito_5_remove012_SCT_MAST%>%group_by(cluster)%>%arrange(desc(avg_log2FC))%>%dplyr::slice(1:100)%>%ungroup()

write_csv(manuscript_marker_genes,file="manuscript_marker_genes.csv")



################################################
# Xist plot
################################################
load("final_Seurat_object.RData")
dat<-final_Seurat_object
male_dat<-subset(dat,sex=="M")

Xist<-male_dat@assays$RNA$data["Xist",]

summary(Xist)

Xist<-male_dat@assays$RNA$data["Xist",]
Ddx3y<-male_dat@assays$RNA$data["Ddx3y",]
ct<-male_dat@meta.data$cell_type

pdf(file="manuscript_Xist_plot.pdf",width=8,height=6)
DefaultAssay(male_dat)<-"RNA"
print(FeaturePlot(male_dat,features="Xist"))
print(DotPlot(male_dat, features = "Xist",assay="RNA") + RotatedAxis())
dev.off()



##############################
# GO Term analysis
##############################
library(Scillus)
library(fgsea)
library(org.Mm.eg.db)
library(clusterProfiler)
library(msigdf)

pdf("GO_plots_to_compare_cts2.pdf"
    ,width=10,height=6)
n_ct<-length(levels(find_all_markers$cluster))
j<-0
for(ct_name in levels(find_all_markers$cluster)){
  j<-j+1
  for(ontology in c("BP","MF","CC")){
    print(plot_cluster_go(find_all_markers,cluster_name=ct_name,org="mouse",ont=ontology)+
            ggtitle(paste0(ct_name," ",ontology," ontology"))+
            xlab("Gene Ontology Term"))
  }
  print(paste0(j," of ",n_ct," done"))
}
dev.off()

ct_order<-levels(Idents(dat))
levs<-paste0(ct_order,"_As")
pdf("GO_plots_As_vs_Control.pdf"
    ,width=10,height=6)
n_ct<-length(levs)
j<-0
for(ct_name in levs){
  if(grepl("_As",ct_name)){
    j<-j+1
    for(ontology in c("BP","MF","CC")){
      print(plot_cluster_go(find_all_markers_As_vs_Control,cluster_name=ct_name,org="mouse",ont=ontology)+
              ggtitle(paste0(ct_name," ",ontology," ontology"))+
              xlab("Gene Ontology Term"))
    }
    print(paste0(j," of ",n_ct," done"))
  }
  
}
dev.off()
