######### 0. Library #####
library(ggpubr)
library(dplyr)
library(Seurat)
library(ggsci)
library(CellChat)

library(GSVA)
library("ggstatsplot")
######### 1. Data import (In server) ########
BCC_sc <- Read10X_h5("/public/home/ac8q5y82z9/data/0.raw/SiT_validation/validation/no_treatment/BCC_GSE139829_expression.h5")  
BCC_sc <- CreateSeuratObject(BCC_sc,project = "BCC") #后面就可以单细胞处理的标准流程啦
BCC_sc_metadata <- read.table("/public/home/ac8q5y82z9/data/0.raw/SiT_validation/validation/no_treatment/BCC_GSE139829_CellMetainfo_table.tsv",sep = "\t",row.names = 1,header = T)

######### 2. Add metadata (In server) #######
BCC_sc <- AddMetaData(
  object = BCC_sc,
  metadata = BCC_sc_metadata
)

######### 3. Scale and PCA (In server) ########
BCC_sc_feature <- FindVariableFeatures(BCC_sc)

BCC_sc_scaled <- ScaleData(BCC_sc_feature)
BCC_sc_PCA <- RunPCA(BCC_sc_scaled,features = VariableFeatures(object = BCC_sc_scaled))

BCC_sc_PCA_umap <- BCC_sc_PCA %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors( dims = 1:10) %>%
  FindClusters( dims = 1:10,resolution = 0.5) %>% 
  identity()
BCC_sc_PCA_tsne <- BCC_sc_PCA_umap %>% 
  RunTSNE( dims = 1:10)

#### save
saveRDS(BCC_sc_PCA_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_sc_PCA_tsne")
BCC_sc_PCA_tsne <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_aPD1_sc_PCA_tsne")
BCC_sc_PCA_tsne <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_aPD1_sc_PCA_tsne")
BCC_sc_PCA_tsne <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/BCC_aPD1_sc_PCA_tsne")


######### 4.1 Extract macrophage and cluster again (In server) ######
mcells_only_sce_anno_ori <- subset(BCC_sc_PCA_tsne , subset = Celltype..major.lineage. == "Mono/Macro")

mcells_only_sce_anno <- SCTransform(mcells_only_sce_anno_ori) %>% RunPCA() 
pc.num=1:10
mcells_only_sce_anno <- RunUMAP(mcells_only_sce_anno, dims=pc.num) %>%
  FindNeighbors( dims=pc.num) %>% 
  FindClusters(resolution=0.8)   
mcells_only_sce_anno@active.assay = "SCT"
saveRDS(mcells_only_sce_anno,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_mcells_only_sce_anno_original.rds")
mcells_only_sce_anno <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_mcells_only_sce_anno.rds")

mcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/BCC_mcells_only_sce_anno.rds")
######### 4.2 Add C1QC label ####
## only choose cluster 8 as C1QC+ TAM
mcells_only_sce_anno$TAM_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(2,3,6,7,10,11), "FCN1+ TAM",
                                         ifelse(mcells_only_sce_anno$seurat_clusters %in% c(8), "C1QC+ TAM", "Other TAM"))

# mcells_only_sce_anno$tam_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(1,4),"C1QC+ TAM",
#                                          ifelse(mcells_only_sce_anno$seurat_clusters %in% c(0,5,6),"FCN1+ TAM",
#                                                 ifelse(mcells_only_sce_anno$seurat_clusters %in% c(1),"Other TAM","deleted")))
# 


######### 5. Add label into original sc #######

TAM_label <- mcells_only_sce_anno@meta.data[,"TAM_label"]
names(TAM_label ) <-  rownames(mcells_only_sce_anno@meta.data)
BCC_sc_TAM_labled <- AddMetaData(
  object = BCC_sc_PCA_tsne,
  metadata = as.data.frame(TAM_label)
)

head(BCC_sc_TAM_labled@meta.data)
BCC_sc_TAM_labled$Celltype..original. <- ifelse(is.na(BCC_sc_TAM_labled$TAM_label), BCC_sc_TAM_labled$Celltype..original., BCC_sc_TAM_labled$TAM_label)
table(BCC_sc_TAM_labled$Celltype..original.)

## In local
saveRDS(BCC_sc_TAM_labled,"~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled.rds" )
BCC_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled.rds")

## In sever
saveRDS(BCC_sc_TAM_labled,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_sc_TAM_labled.rds" )
BCC_sc_TAM_labled <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_sc_TAM_labled.rds" )
######### 6.1 Marker #########
## test c1qc and fcn1 diff markers
c1qc_TAM_markers <- FindMarkers(BCC_sc_TAM_labled, group.by = "Celltype..original.",ident.1 = "C1QC+ TAM", ident.2 = "FCN1+ TAM")
c1qc_TAM_markers$gene <- rownames(c1qc_TAM_markers)
c1qc_TAM_markers = c1qc_TAM_markers%>% select(gene, everything()) %>% subset(p_val<0.05)
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
top50 = top_n(c1qc_TAM_markers ,n = 50, wt = abs(avg_log2FC))
write.csv(c1qc_TAM_markers, "~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/c1qc_TAM_markers.csv", row.names = F)
write.csv(top50, "~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/top50_c1qc_TAM_markers.csv", row.names = F)

## test c1qc diff other all cluster markers

######### 6.2 Marker visulization #########
## 1) BCC_sc_PCA_tsne
pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_marker_plot.pdf",width = 10)
DimPlot(BCC_sc_PCA_tsne, reduction = "tsne", group.by = "Celltype..original.")
FeaturePlot(BCC_sc_PCA_tsne, reduction = "tsne",features = c("CXCL13", "C1QC", "FCN1"))
VlnPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Celltype..major.lineage.")
VlnPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Celltype..minor.lineage.")
VlnPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Celltype..original.")
VlnPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Response" )
DotPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),group.by = "Celltype..major.lineage.")
DotPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),group.by = "Response")
DotPlot(BCC_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),group.by = "Celltype..original.")
dev.off()

### 2）mcells_only_sce_anno
pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_mcell_marker_plot_original.pdf",width = 8)
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "seurat_clusters", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Celltype..major.lineage.", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Celltype..minor.lineage.", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Response", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "TimePoint", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Celltype..original.", label = F) 
#DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "TAM_label", label = F) 

FeaturePlot(mcells_only_sce_anno, reduction = "tsne",features = c("FCN1", "C1QC", "CD74","CCL5"))
VlnPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),pt.size = 0,group.by ="seurat_clusters")
VlnPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),pt.size = 0,group.by ="Response")
#VlnPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),pt.size = 0,group.by ="Celltype..minor.lineage.")
#VlnPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),pt.size = 0,group.by ="TAM_label")
DotPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),group.by = "seurat_clusters")
#DotPlot(mcells_only_sce_anno,features =  c("FCN1", "C1QC", "CD74"),group.by = "TAM_label")

dev.off()

### 3）BCC_sc_TAM_labled
#pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/BCC_TAM_labled_marker_plot.pdf",width = 10)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/original_BCC_TAM_labled_marker_plot.pdf",width = 10)

DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "Celltype..major.lineage.", label = F) 
DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "Celltype..minor.lineage.", label = F) 
DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "Response", label = F) 
DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "TimePoint", label = F) 
DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "Celltype..original.", label = F) 

FeaturePlot(BCC_sc_TAM_labled, reduction = "tsne",features = c("FCN1", "C1QC", "CD74","CXCL13"))
VlnPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),pt.size = 0,group.by ="TimePoint")
VlnPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),pt.size = 0,group.by ="Response")
VlnPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),pt.size = 0,group.by ="Celltype..major.lineage.")

DotPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),group.by = "Celltype..major.lineage.")
DotPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),group.by = "Celltype..minor.lineage.")
DotPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),group.by = "Celltype..original.")
DotPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),group.by = "Response")
DotPlot(BCC_sc_TAM_labled,features =  c("FCN1", "C1QC", "CXCL13"),group.by = "TimePoint")

dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/original_BCC_marker_dotplot.pdf",width = 8)
DimPlot(BCC_sc_TAM_labled, reduction = "tsne", group.by = "Celltype..original.", label = F) 

DotPlot(BCC_sc_TAM_labled,features = c("FCN1", "C1QC", "CXCL13"), group.by = "Celltype..original.") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
FeaturePlot(BCC_sc_TAM_labled,reduction = "tsne",features = c("FCN1", "C1QC", "CD74","CXCL13"))

dev.off()
#### conclusion:  original better
BCC_sc_TAM_labled@meta.data$Celltype..original. <- gsub("Macrophages","DCs",BCC_sc_TAM_labled@meta.data$Celltype..original.)

saveRDS(BCC_sc_TAM_labled,"~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled.rds" )
BCC_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled.rds")
######### 7. Ratio ######

### (1) boxplot （pre R vs pre NR vs post R vs post NR) ####
group_response_compare <- function(metadata, celltype){
  boxplot_df_all <- prop.table(table(metadata$Celltype..original., metadata$Sample), margin = 2)
  boxplot_df_all <- boxplot_df_all[celltype,]
  boxplot_df_all <- data.frame("celltype" = boxplot_df_all, "Sample" = names(boxplot_df_all))
  boxplot_df_all$Response <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",3)
  boxplot_df_all$Patient <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",1)
  boxplot_df_all$Treatment <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",2)
  boxplot_df_all$celltype <- as.numeric(boxplot_df_all$celltype)
  
  ## response vs patient
  response_patient_df <-boxplot_df_all[,c("Response", "Patient")]
  response_patient_df <- response_patient_df[response_patient_df$Response != "None",]
  boxplot_df_all <- merge(boxplot_df_all,response_patient_df,by = "Patient")
  boxplot_df_all$Group <- paste0(boxplot_df_all$Treatment,"_",boxplot_df_all$Response.y)
  
  ggboxplot(boxplot_df_all ,x = "Group", y = "celltype" ,fill = "Group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
    stat_compare_means(comparisons = list(c("pre_R", "pre_NR"),
                                                            c("pre_R", "post_R"),
                                                            c("pre_NR", "post_NR"))) +
    ylab(celltype)
}
BCC_sc_TAM_labled@meta.data$Sample <- paste0(BCC_sc_TAM_labled$Patient,"_", BCC_sc_TAM_labled$TimePoint,"_", BCC_sc_TAM_labled$Response)

### ratio in all cells
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/w_test_response_comparison_ratio_boxplot.pdf",width = 5, height = 5)
group_response_compare(BCC_sc_TAM_labled, "C1QC+ TAM")
group_response_compare(BCC_sc_TAM_labled, "CD8_ex_T_cells")
group_response_compare(BCC_sc_TAM_labled, "FCN1+ TAM")
group_response_compare(BCC_sc_TAM_labled, "CD4_T_cells")
group_response_compare(BCC_sc_TAM_labled, "CD8_mem_T_cells")
group_response_compare(BCC_sc_TAM_labled, "CD8_act_T_cells")
group_response_compare(BCC_sc_TAM_labled, "Other TAM")
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/Celltype_malignancy_response_comparison_ratio_boxplot.pdf",width = 5, height = 5)

group_response_compare(BCC_sc_TAM_labled, "Immune cells" )
group_response_compare(BCC_sc_TAM_labled, "Malignant cells")

dev.off()


### ratio in imm cells
BCC_sc_TAM_labled_imm <- BCC_sc_TAM_labled@meta.data[BCC_sc_TAM_labled$Celltype..malignancy. == "Immune cells",]
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/response_comparison_ratio_imm_cells_boxplot.pdf",width = 5, height = 5)
group_response_compare(BCC_sc_TAM_labled_imm, "C1QC+ TAM")
group_response_compare(BCC_sc_TAM_labled_imm, "CD8_ex_T_cells")
group_response_compare(BCC_sc_TAM_labled_imm, "FCN1+ TAM")
group_response_compare(BCC_sc_TAM_labled_imm, "CD4_T_cells")
group_response_compare(BCC_sc_TAM_labled_imm, "CD8_mem_T_cells")
group_response_compare(BCC_sc_TAM_labled_imm, "CD8_act_T_cells")
dev.off()


### (2) BCC ggbarstats #####
##  only TAM
BCC_sc_TAM_labled_TAM <-  BCC_sc_TAM_labled@meta.data[BCC_sc_TAM_labled$Celltype..original. %in% c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_sc_TAM_labled_TAM_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(BCC_sc_TAM_labled_TAM, Stage, Celltype..original., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(BCC_sc_TAM_labled_TAM,  Celltype..original., Stage,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()

## CD8+ T cells
BCC_sc_TAM_labled_CD8 <-  BCC_sc_TAM_labled@meta.data[BCC_sc_TAM_labled$Celltype..original. %in% c("CD4Tconv", "CD8T", ' CD8_ex_T_cells'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_sc_TAM_labled_CD8_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(BCC_sc_TAM_labled_CD8, Stage, Celltype..original., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(BCC_sc_TAM_labled_CD8,  Celltype..original., Stage,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()


# group_response_compare_ <- function(metadata, celltype){
#   boxplot_df_all <- prop.table(table(metadata$seurat_clusters, metadata$Sample), margin = 2)
#   boxplot_df_all <- boxplot_df_all[celltype,]
#   boxplot_df_all <- data.frame("celltype" = boxplot_df_all, "Sample" = names(boxplot_df_all))
#   boxplot_df_all$Response <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",3)
#   boxplot_df_all$Patient <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",1)
#   boxplot_df_all$Treatment <- sapply(strsplit(boxplot_df_all$Sample,"_"),"[[",2)
#   boxplot_df_all$celltype <- as.numeric(boxplot_df_all$celltype)
#   
#   ## response vs patient
#   response_patient_df <-boxplot_df_all[,c("Response", "Patient")]
#   response_patient_df <- response_patient_df[response_patient_df$Response != "None",]
#   boxplot_df_all <- merge(boxplot_df_all,response_patient_df,by = "Patient")
#   boxplot_df_all$Group <- paste0(boxplot_df_all$Treatment,"_",boxplot_df_all$Response.y)
#   
#   ggboxplot(boxplot_df_all ,x = "Group", y = "celltype" ,fill = "Group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
#     stat_compare_means(method = "t.test",comparisons = list(c("pre_R", "pre_NR"),
#                                                             c("pre_R", "post_R"),
#                                                             c("pre_NR", "post_NR"))) +
#     ylab(celltype)
# }
# 
# pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_mcells_response_compare.pdf")
# group_response_compare_(mcells_only_sce_anno, "1")
# group_response_compare_(mcells_only_sce_anno, "2")
# group_response_compare_(mcells_only_sce_anno, "3")
# group_response_compare_(mcells_only_sce_anno, "4")
# group_response_compare_(mcells_only_sce_anno, "5")
# group_response_compare_(mcells_only_sce_anno, "6")
# group_response_compare_(mcells_only_sce_anno, "7")
# group_response_compare_(mcells_only_sce_anno, "8")
# group_response_compare_(mcells_only_sce_anno, "9")
# group_response_compare_(mcells_only_sce_anno, "10")
# group_response_compare_(mcells_only_sce_anno, "11")
# dev.off()





###（3）CD74 exp vs CD8 effector ratio vs C1QC TAM ratio ######
cell_ratio_sample_total <- prop.table(table(BCC_sc_TAM_labled$Celltype..original.,BCC_sc_TAM_labled$Sample), margin = 2)

### only C1QC TAM
c1qc_only <- subset(BCC_sc_TAM_labled, subset = Celltype..original. == "C1QC+ TAM")

expr_c1qc_only <- AverageExpression(c1qc_only, assays = "RNA", slot = "data",group.by = "Sample")
expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]
expr_c1qc_only <- data.frame("CD74_expression" = expr_c1qc_only, "Patient" = names(expr_c1qc_only))

dotplot_5 <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                        "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                       sample = colnames(cell_ratio_sample_total))

## add info
dotplot_5$Patient <- sapply(strsplit(dotplot_5$sample,"_"),"[[",1)
dotplot_5$Response <- ifelse(dotplot_5$Patient %in% pre_R_patient ,"R","NR")
dotplot_5$Treatment <- sapply(strsplit(dotplot_5$sample,"_"),"[[",2)
dotplot_5$Group <- paste0(dotplot_5$Treatment,"_",dotplot_5$Response)

## merge cd74 expression
dotplot_5 <- merge(dotplot_5, expr_c1qc_only, by.x = "sample" , by.y = "Patient",all.x = T)
dotplot_5$CD74_expression <- ifelse(is.na(dotplot_5$CD74_expression), 0,dotplot_5$CD74_expression)

p5 <- ggplot(dotplot_5, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=CD74_expression, color=Group)) +
  geom_point() + scale_size(range = c(1, 10)) + scale_color_viridis(discrete = T) + theme_classic() + theme(legend.position = "top")


pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_ratio_cd74_expression_point_plot.pdf",width= 8) 
p5
dev.off()
### (4) CD74 expression boxplot ######

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/w_test_BCC_cd74_expression_boxplot.pdf")
ggboxplot(dotplot_5 ,x = "Group", y = "CD74_expression" ,fill = "Group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("pre_R", "pre_NR"),
                                                          c("pre_R", "post_R"),
                                                          c("pre_NR", "post_NR")) )
dev.off()


### (5) All celltype vs patient ####
cell_ratio_sample <- prop.table(table(BCC_sc_TAM_labled$Celltype..original., BCC_sc_TAM_labled$Patient), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled_sample_all_ratio.pdf",width = 5, height = 5)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() + scale_fill_viridis(option="b", discrete = T) +
  labs(x='Sample',y = 'Ratio') + coord_flip() +
   theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "bottom")  
dev.off()

### (7) All celltype vs response ####
BCC_sc_TAM_labled$Group <- ifelse(BCC_sc_TAM_labled$Treatment == "None" & BCC_sc_TAM_labled$Patient %in% pre_NR_patient,"Pre_NR",
                                  ifelse(BCC_sc_TAM_labled$Treatment == "None" & BCC_sc_TAM_labled$Patient %in% pre_R_patient,"Pre_R",
                                         ifelse(BCC_sc_TAM_labled$Response == "NR","Podt_NR",
                                                ifelse(BCC_sc_TAM_labled$Response == "R","Podt_R","None"))))
                                                
table(BCC_sc_TAM_labled$Group)       
cell_ratio_sample <- prop.table(table(BCC_sc_TAM_labled$Celltype..original., BCC_sc_TAM_labled$Group), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)


pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled_group.pdf",width = 8, height = 5)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() + scale_fill_viridis(option="A", discrete = T,direction = -1) +
  labs(x='Sample',y = 'Ratio') +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))  
dev.off()

### (7) Response vs All celltype  ####
table(BCC_sc_TAM_labled$Group)       
cell_ratio_sample <- prop.table(table( BCC_sc_TAM_labled$Group, BCC_sc_TAM_labled$Celltype..original.), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)


pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_sc_TAM_labled_group_cell.pdf",width = 5, height = 4)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() + scale_fill_viridis(option="A", discrete = T,direction = -1) +
  labs(x='Sample',y = 'Ratio') + coord_flip() +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "bottom")  
dev.off()
### (8) Celltype piechart ######
library(lessR)
BCC_sc_TAM_labled_df <- as.data.frame(BCC_sc_TAM_labled@meta.data)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/celltype_ratio_piechart.pdf")
PieChart(Celltype..original., data = BCC_sc_TAM_labled_df,
         main = NULL, fill = "blues") 
dev.off()
######### 9. BCC Cellchat all ######
library(CellChat)
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/cellchat_all/")

## 10.1 创建cellchat对象
cellchat <- createCellChat(BCC_sc_TAM_labled , group.by = "Celltype..original.")
cellchat <- setIdent(cellchat, ident.use = "Celltype..original.") # set "labels" as default cell identity

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
## 10.2 设置参考数据库
cellchat@DB <- CellChatDB.human
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part2: Inference of cell-cell conmmunication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "BCC_cellchat.rds")
#cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/BCC_cellchat.rds")
#cellchat <- readRDS("BCC_cellchat.rds")
#### 10.3 可视化 
## (0) save csv
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv")

df.netp <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")
## (1) netVisual_circle
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("netVisual_circle.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
pdf("netVisual_circle_sep.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

### (2) netVisual_heatmap

par(mfrow=c(1,1))
pdf("netVisual_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = "APP", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "CXCL", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "CCL", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "IL1", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "MHC-I", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "MHC-II", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "VEGF", color.heatmap = "Reds")

netVisual_heatmap(cellchat, color.heatmap = "Reds")
dev.off()

pdf("weight_netVisual_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = "APP", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "CXCL", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "CCL", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "IL1", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "MHC-I", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "MHC-II", color.heatmap = "Reds",measure = "weight")
netVisual_heatmap(cellchat, signaling = "VEGF", color.heatmap = "Reds",measure = "weight")

netVisual_heatmap(cellchat, color.heatmap = "Reds",measure = "weight")
dev.off()

####
pathways.show <- c("CXCL", "CCL","MIF","MHC-I","MHC-II","APP") 
pdf("sig_pathway_L-R.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("MIF_sig_pathway_L-R.pdf")
netAnalysis_contribution(cellchat, signaling = "MIF")
dev.off()

### (4) netVisual_embedding 信号网络聚类
# 按功能相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
p = netVisual_embedding(cellchat, type = "functional")
ggsave("custer_pathway_function.pdf", p, width =5, height = 4)
p = netVisual_embeddingZoomIn(cellchat, type = "functional")
ggsave("custer_pathway_function2.pdf", p, width = 8, height = 6)

# 按结构相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
p = netVisual_embedding(cellchat, type = "structural")
ggsave("custer_pathway_structure.pdf", p, width = 9, height = 6)
p = netVisual_embeddingZoomIn(cellchat, type = "structural")
ggsave("custer_pathway_structure2.pdf", p, width = 8, height = 6)

##### only labelede cd8 cxclq3 and c1qc tam
p = netVisual_embedding(cellchat, type = "structural",pathway.labeled = c1qc_cd8_scores$pathway_name)
ggsave("labeled_custer_pathway_structure2.pdf", p, width = 5, height = 5)

p = netVisual_embedding(cellchat, type = "functional",pathway.labeled = c1qc_cd8_scores$pathway_name)
ggsave("labeled_custer_pathway_functional.pdf", p, width = 5, height = 5)


#### (5) 单个信号通路
cellchat@netP$pathways
pdf("sep_netVisual_aggregate_chord.pdf")
netVisual_aggregate(cellchat,signaling = "MIF", layout = "chord") 
netVisual_aggregate(cellchat,signaling = "MIF", layout = "chord",
                    sources.use = c(4,5), targets.use = c(2,7,9)) 
netVisual_aggregate(cellchat,signaling = "MIF", layout = "chord",
                    sources.use = c(5), targets.use = c(2,7,9)) 

dev.off()
#### (6) 单 配体受体
pairLG.MIF <- extractEnrichedLR(cellchat, signaling = "MIF")
pdf("MIF_netVisual_aggregate_chord.pdf")
netVisual_individual(cellchat,signaling = "MIF",pairLR.use = pairLG.MIF$interaction_name[1], 
                     layout = "chord", sources.use = c(5), targets.use = c(2,7,9))
netVisual_individual(cellchat,signaling = "MIF",pairLR.use = pairLG.MIF$interaction_name[2], 
                     layout = "chord", sources.use = c(5), targets.use = c(2,7,9))

netVisual_individual(cellchat,signaling = "MIF",pairLR.use = pairLG.MIF$interaction_name[1],
                     sources.use = c(5), targets.use = c(2,7,9),layout = "chord")
netVisual_individual(cellchat,signaling = "MIF",pairLR.use = pairLG.MIF$interaction_name[2],
                     sources.use = c(5), targets.use = c(2,7,9),layout = "chord")

dev.off()



####### (7) 多个配体受体介导的细胞互作关系可视化
## cd8 effector to all TAM
p1 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(2,7,9),
                       remove.isolate = F)
p2 <- netVisual_bubble(cellchat,  sources.use = c(5), targets.use = c(2,7,9),
                       remove.isolate = F, signaling = "MIF")



pdf("all_MIF_cd8_effector_all_TAM.pdf")
p1
p2
dev.off()

##### Lolipop
c1qc_cd8_scores <- subsetCommunication(cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8_ex_T_cells" & c1qc_cd8_scores$target == "C1QC+ TAM",]

pdf("c1qc_cd8_interaction.pdf",width = 10)
ggdotchart(c1qc_cd8_scores, y = "prob", x ="interaction_name_2", color = "interaction_name_2",dot.size = 6,
           add = "segments") + scale_color_viridis(discrete = T) + theme(legend.position = "none",axis.text.x = element_text(angle = 35)) 

dev.off()


####### 9. Overlap pathway #####
### BCC
BCC_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/cellchat_all/cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(BCC_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8_ex_T_cells" & c1qc_cd8_scores$target == "C1QC+ TAM",]
BCC_pathway <- unique(c1qc_cd8_scores$pathway_name)

### UVM
UVM_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/cellchat_all/UVM_cellchat.rds")

c1qc_cd8_scores <- subsetCommunication(UVM_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ CXCL13+ Tex" & c1qc_cd8_scores$target == "C1QC+ TAM",]
UVM_pathway <- unique(c1qc_cd8_scores$pathway_name)

### SiT
sit_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(sit_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ Effctor T" & c1qc_cd8_scores$target == "C1QC+ TAM",]
sit_pathway <- unique(c1qc_cd8_scores$pathway_name)

list_pathway <- list("BCC" =BCC_pathway, "UVM" = UVM_pathway, "ESCC" = sit_pathway)

venn(list_pathway)
venn.plot <- venn.diagram(
  x = list_pathway,
  filename = "~/onedrive/Work/phD/phd_project/SiT/results/validation/venn.PDF",
  col = "black",
  #lty = "dotted",
  fill = c("#300A5BFF", "#C43C4EFF", "#FCA007FF"),
  alpha = 0.9,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,#内标签的字体大小
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "outer",#设置标签在圆外面
  #cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2,#外标签的字体大小
  cat.fontfamily = "serif",
  cat.dist = c(0.05, 0.05, 0.05),#相对圆圈的位置
  cat.pos = c(-20,20,180)  #相对12点方向旋转的角度
)

intersect(intersect(BCC_pathway,UVM_pathway),sit_pathway)
setdiff(sit_pathway,BCC_pathway)

setdiff(intersect(sit_pathway,BCC_pathway),sit_pathway)
######### 9. BCC Sep Cellchat only before (archive) ########
response_patient_df <- BCC_sc_TAM_labled@meta.data[,c("Response", "Patient")]
response_patient_df <- response_patient_df[response_patient_df$Response != "None",]
response_patient_df <- response_patient_df[!duplicated(response_patient_df),]
pre_R_patient <- response_patient_df[response_patient_df$Response == "R", "Patient"]
pre_NR_patient <- response_patient_df[response_patient_df$Response == "NR", "Patient"]

###  extract only pre treatment 
pre_BCC_sc_TAM_labled <- subset(BCC_sc_TAM_labled, subset = TimePoint == "pre")

R_BCC_sce <- subset(pre_BCC_sc_TAM_labled, subset = Patient %in% pre_R_patient)
NR_BCC_sce <- subset(pre_BCC_sc_TAM_labled, subset = Patient %in% pre_NR_patient)

setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_before_sep")

### 11.1 create cellchat 
R_cellchat <- createCellChat(R_BCC_sce, group.by = "Celltype..original.")
NR_cellchat <- createCellChat(NR_BCC_sce, group.by = "Celltype..original.")

### 11.2 构建细胞通讯网络
## R
R_cellchat@DB <- CellChatDB.human
R_cellchat <- subsetData(R_cellchat) # This step is necessary even if using the whole database
R_cellchat <- identifyOverExpressedGenes(R_cellchat)
R_cellchat <- identifyOverExpressedInteractions(R_cellchat)

R_cellchat <- computeCommunProb(R_cellchat)
R_cellchat <- filterCommunication(R_cellchat, min.cells = 10)
R_cellchat <- computeCommunProbPathway(R_cellchat)
R_cellchat <- aggregateNet(R_cellchat)

saveRDS(R_cellchat, "R_cellchat.rds")
R_cellchat <- readRDS("R_cellchat.rds")

## NR
NR_cellchat@DB <- CellChatDB.human
NR_cellchat <- subsetData(NR_cellchat) # This step is necessary even if using the whole database
NR_cellchat <- identifyOverExpressedGenes(NR_cellchat)
NR_cellchat <- identifyOverExpressedInteractions(NR_cellchat)

NR_cellchat <- computeCommunProb(NR_cellchat)
NR_cellchat <- filterCommunication(NR_cellchat, min.cells = 10)
NR_cellchat <- computeCommunProbPathway(NR_cellchat)
NR_cellchat <- aggregateNet(NR_cellchat)

saveRDS(NR_cellchat, "NR_cellchat.rds")
NR_cellchat <- readRDS("NR_cellchat.rds")

### 11.3 Merge CellChat object
R_NR_cellchat_list <- list(R = R_cellchat, NR = NR_cellchat)
R_NR_merged_cellchat <- mergeCellChat(R_NR_cellchat_list, add.names = names(R_NR_cellchat_list),
                                      cell.prefix = T)


saveRDS(R_NR_merged_cellchat,"R_NR_merged_cellchat")

### 11.4 可视化

# （1） 通讯数量与强度对比
#gg1 <- compareInteractions(R_NR_merged_cellchat, show.legend = F , group = c(1,2))
gg2 <- compareInteractions(R_NR_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("R_NR_compareInteractions.pdf")
gg1 + gg2
dev.off()

#gg1 <- compareInteractions(R_NR_merged_cellchat, show.legend = F,group = c(1,2))
gg2 <- compareInteractions(R_NR_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("R_NR_compareInteractions.pdf")
gg1 + gg2
dev.off()


# (2) 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_groups_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(R_NR_merged_cellchat, weight.scale = T, measure = "weight")

dev.off()

pdf("TAM_CD8_tex_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(R_NR_merged_cellchat, weight.scale = T, measure = "weight",sources.use = c(7), targets.use = c(3,11,16))

dev.off()

# (3) Heatmap
pdf("R_NR_netVisual_heatmap.pdf")
netVisual_heatmap(R_NR_merged_cellchat, measure = "weight")
dev.off()


# (4) 特定信号通路的对比
weight.max <- getMaxWeight(R_NR_cellchat_list , attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_netVisual_circle.pdf")
for (i in 1:length(R_NR_cellchat_list )) {
  netVisual_circle(R_NR_cellchat_list [[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", 
                                                            names(R_NR_cellchat_list )[i]))
}
dev.off()

pdf("MIF_netVisual_circle.pdf")
for (i in 1:length(R_NR_cellchat_list )) {
  netVisual_aggregate(R_NR_cellchat_list[[i]], layout = "circle", signaling = c("MIF"),
                      label.edge= F, edge.weight.max = weight.max[2], 
                      edge.width.max = 12, signaling.name = paste0("MIF  - ", 
                                                                   names(R_NR_cellchat_list )[i]))
}
dev.off()

# (5) 通路信号强度对比分析

pdf("R_NR_rankNet.pdf",height = 10)
rankNet(R_NR_merged_cellchat, mode= "comparison",stacked = T)
dev.off()


# (6) 配体-受体对比分析
pdf("R_NR_TAM_cd8_total_netVisual_bubble.pdf")
netVisual_bubble(R_NR_merged_cellchat, sources.use = c(6), targets.use = c(2,3,4),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in R", max.dataset = 1)
netVisual_bubble(R_NR_merged_cellchat, sources.use = c(6), targets.use = c(2,3,4),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in NR", max.dataset = 2)
dev.off()













######### 10. BCC Sep Cellchat before vs after R (archive) ########
response_patient_df <- BCC_sc_TAM_labled@meta.data[,c("Response", "Patient")]
response_patient_df <- response_patient_df[response_patient_df$Response != "None",]
response_patient_df <- response_patient_df[!duplicated(response_patient_df),]
pre_R_patient <- response_patient_df[response_patient_df$Response == "R", "Patient"]
pre_NR_patient <- response_patient_df[response_patient_df$Response == "NR", "Patient"]

###  extract only pre treatment 
pre_BCC_sc_TAM_labled <- subset(BCC_sc_TAM_labled, subset = TimePoint == "pre")
post_BCC_sc_TAM_labled <- subset(BCC_sc_TAM_labled, subset = TimePoint == "post")

pre_R_sce <- subset(pre_BCC_sc_TAM_labled, subset = Patient %in% pre_R_patient)
pre_NR_sce <- subset(pre_BCC_sc_TAM_labled, subset = Patient %in% pre_NR_patient)

post_R_sce <- subset(post_BCC_sc_TAM_labled, subset = Patient %in% pre_R_patient)
post_NR_sce <- subset(post_BCC_sc_TAM_labled, subset = Patient %in% pre_NR_patient)

setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_pre_post")

### 11.1 create cellchat 
cellchat_pipeline <- function(seurat,groupby, rds_name){
  cellchat <- createCellChat(seurat, group.by = groupby)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, rds_name)
  return(cellchat)
}

pre_R_cellchat <- cellchat_pipeline(pre_R_sce, "Celltype..original.","pre_R_cellchat.rds")
pre_NR_cellchat <- cellchat_pipeline(pre_NR_sce, "Celltype..original.","pre_NR_cellchat.rds")
post_R_cellchat <- cellchat_pipeline(post_R_sce, "Celltype..original.","post_R_cellchat.rds")
post_NR_cellchat <- cellchat_pipeline(post_NR_sce, "Celltype..original.","post_NR_cellchat.rds")


### 11.3 Merge CellChat object
R_pre_post_cellchat_list <- list(pre = pre_R_cellchat, post = post_R_cellchat)
R_pre_post_merged_cellchat <- mergeCellChat(R_pre_post_cellchat_list, add.names = names(R_pre_post_cellchat_list),
                                            cell.prefix = T)

NR_pre_post_cellchat_list <- list(pre = pre_NR_cellchat, post = post_NR_cellchat)
NR_pre_post_merged_cellchat <- mergeCellChat(NR_pre_post_cellchat_list, add.names = names(NR_pre_post_cellchat_list),
                                             cell.prefix = T)

saveRDS(R_pre_post_merged_cellchat,"R_pre_post_merged_cellchat.rds")
saveRDS(NR_pre_post_merged_cellchat,"NR_pre_post_merged_cellchat.rds")

### 11.4 可视化

# （1） 通讯数量与强度对比
gg1 <- compareInteractions(R_pre_post_merged_cellchat, show.legend = F , group = c(1,2), measure = "weight")
gg2 <- compareInteractions(NR_pre_post_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("pre_post_compareInteractions.pdf")
gg1 + gg2
dev.off()


# (2) 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf("pre_post_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(R_pre_post_merged_cellchat, weight.scale = T, measure = "weight")
netVisual_diffInteraction(NR_pre_post_merged_cellchat, weight.scale = T, measure = "weight")

dev.off()

pdf("TAM_CD8_tex_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(R_pre_post_merged_cellchat, weight.scale = T, measure = "weight",sources.use = c(7), targets.use = c(3,11,16))
netVisual_diffInteraction(NR_pre_post_merged_cellchat, weight.scale = T, measure = "weight",sources.use = c(7), targets.use = c(3,11,16))

dev.off()

# (3) Heatmap
pdf("pre_post_netVisual_heatmap.pdf")
netVisual_heatmap(R_pre_post_merged_cellchat, measure = "weight")
netVisual_heatmap(NR_pre_post_merged_cellchat, measure = "weight")

dev.off()


# (4) 特定信号通路的对比
weight.max <- getMaxWeight(R_pre_post_merged_cellchat , attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("pre_post_R_netVisual_circle.pdf")
for (i in 1:length(R_pre_post_merged_cellchat )) {
  netVisual_circle(R_pre_post_merged_cellchat [[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", 
                                                            names(R_pre_post_merged_cellchat )[i]))
}
dev.off()

pdf("MIF_netVisual_circle.pdf")
for (i in 1:length(R_pre_post_merged_cellchat )) {
  netVisual_aggregate(R_pre_post_merged_cellchat[[i]], layout = "circle", signaling = c("MIF"),
                      label.edge= F, edge.weight.max = weight.max[2], 
                      edge.width.max = 12, signaling.name = paste0("MIF  - ", 
                                                                   names(R_pre_post_merged_cellchat )[i]))
}
dev.off()

# (5) 通路信号强度对比分析

pdf("R_pre_post_rankNet.pdf",height = 10)
rankNet(R_pre_post_merged_cellchat, mode= "comparison",stacked = T)
dev.off()

pdf("NR_pre_post_rankNet.pdf",height = 10)
rankNet(NR_pre_post_merged_cellchat, mode= "comparison",stacked = T)
dev.off()


# (6) 配体-受体对比分析
pdf("R_pre_post_TAM_cd8_total_netVisual_bubble.pdf")
netVisual_bubble(R_pre_post_merged_cellchat, sources.use = c(7),targets.use = c(3,11),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Pre", max.dataset = 1)
netVisual_bubble(R_pre_post_merged_cellchat, sources.use = c(7),targets.use = c(3,11),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Post", max.dataset = 2)
dev.off()

pdf("NR_pre_post_TAM_cd8_total_netVisual_bubble.pdf")
netVisual_bubble(NR_pre_post_merged_cellchat, sources.use = c(7),targets.use = c(3,11),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Pre", max.dataset = 1)
netVisual_bubble(NR_pre_post_merged_cellchat, sources.use = c(7),targets.use = c(3,11),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Post", max.dataset = 2)
dev.off()

######### 11. BCC pre R vs NR MIF-CD74和 TAM  CD8 CXCL13 ratio correlation (archive)  ########
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_before_sep/")
#BCC_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/BCC_sc_TAM_labled.rds")

#post_NR_cellchat  <- readRDS("post_NR_cellchat.rds")
#pre_NR_cellchat <- readRDS("pre_NR_cellchat.rds")

### pre_R
pre_R_mif_cd74_scores <- subsetCommunication(R_cellchat,slot.name = "netP")
pre_R_mif_cd74_scores <- pre_R_mif_cd74_scores[pre_R_mif_cd74_scores$source == "CD8_ex_T_cells" & pre_R_mif_cd74_scores$target == "C1QC+ TAM",]
pre_R_mif_cd74_scores <- pre_R_mif_cd74_scores[pre_R_mif_cd74_scores$pathway_name == "MIF","prob"]
pre_R_mif_cd74_scores <- 0
### pre_NR
pre_NR_mif_cd74_scores <- subsetCommunication(NR_cellchat,slot.name = "netP")
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$source == "CD8_ex_T_cells" & pre_NR_mif_cd74_scores$target == "C1QC+ TAM",]
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$pathway_name == "MIF","prob"]

## 1) total ratio
## only keep none
BCC_sc_TAM_labled_none <- subset(BCC_sc_TAM_labled, subset = Response == "None")
BCC_sc_TAM_labled_none$Response <- ifelse(BCC_sc_TAM_labled_none$Patient %in% pre_NR_patient, "NR", "R")
cell_ratio_sample_total <- prop.table(table( BCC_sc_TAM_labled_none$Celltype..original. ,BCC_sc_TAM_labled_none$Response), margin = 2)


dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                         cd74_prob = c(pre_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_pre_NR_R_total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


## 2) C1QC / FCN1 TAM ratio
tam_only_metadata <- BCC_sc_TAM_labled_none@meta.data[BCC_sc_TAM_labled_none@meta.data$Celltype..original. %in% c("FCN1+ TAM","C1QC+ TAM"),]
cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..original.,tam_only_metadata$Response), margin = 2)

dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
                             "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                             cd74_prob = c(pre_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_pre_NR_R_FCN1_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)

dev.off()


## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8_ex_T_cells",],
                               cd74_prob = c(pre_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_pre_NR_R_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()

######### 12. BCC pre post R MIF-CD74和 TAM  CD8 CXCL13 ratio correlation (archive) ########
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/BCC_pre_post/")
BCC_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC/BCC_sc_TAM_labled.rds")

post_R_cellchat  <- readRDS("post_R_cellchat.rds")
pre_R_cellchat <- readRDS("pre_R_cellchat.rds")

### post_R
post_R_mif_cd74_scores <- subsetCommunication(post_R_cellchat,slot.name = "netP")
post_R_mif_cd74_scores <- post_R_mif_cd74_scores[post_R_mif_cd74_scores$source == "CD8_ex_T_cells" & post_R_mif_cd74_scores$target == "C1QC+ TAM",]
post_R_mif_cd74_scores <- post_R_mif_cd74_scores[post_R_mif_cd74_scores$pathway_name == "MIF","prob"]

### pre_R
pre_R_mif_cd74_scores <- subsetCommunication(pre_R_cellchat,slot.name = "netP")
pre_R_mif_cd74_scores <- pre_R_mif_cd74_scores[pre_R_mif_cd74_scores$source == "CD8_ex_T_cells" & pre_R_mif_cd74_scores$target == "C1QC+ TAM",]
pre_R_mif_cd74_scores <- pre_R_mif_cd74_scores[pre_R_mif_cd74_scores$pathway_name == "MIF","prob"]
pre_R_mif_cd74_scores = 0
## 1) total ratio
BCC_onlyR <- subset(BCC_sc_TAM_labled , subset = Patient %in% pre_R_patient )
cell_ratio_sample_total <- prop.table(table( BCC_onlyR $Celltype..original. ,BCC_onlyR $TimePoint), margin = 2)

dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                         cd74_prob = c(post_R_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_R_total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


# ## 2) C1QC / FCN1 TAM ratio
# tam_only_metadata <- BCC_sc_TAM_labled@meta.data[BCC_sc_TAM_labled@meta.data$Celltype..original. %in% c("FCN1+ TAM","C1QC+ TAM"),]
# cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..original.,tam_only_metadata$Treatment), margin = 2)
# 
# dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
#                              "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
#                              cd74_prob = c(post_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
# pdf("BCC_NR_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
# ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
#   geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
# dev.off()


## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8_ex_T_cells",],
                               cd74_prob = c(post_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_R_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()



######### 13. BCC pre post NR MIF-CD74和 TAM  CD8 CXCL13 ratio correlation (archive) ########
### post_NR
post_NR_mif_cd74_scores <- subsetCommunication(post_NR_cellchat,slot.name = "netP")
post_NR_mif_cd74_scores <- post_NR_mif_cd74_scores[post_NR_mif_cd74_scores$source == "CD8_ex_T_cells" & post_NR_mif_cd74_scores$target == "C1QC+ TAM",]
post_NR_mif_cd74_scores <- post_NR_mif_cd74_scores[post_NR_mif_cd74_scores$pathway_name == "MIF","prob"]

### pre_NR
pre_NR_mif_cd74_scores <- subsetCommunication(pre_NR_cellchat,slot.name = "netP")
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$source == "CD8_ex_T_cells" & pre_NR_mif_cd74_scores$target == "C1QC+ TAM",]
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$pathway_name == "MIF","prob"]

## 1) total ratio
BCC_onlyNR <- subset(BCC_sc_TAM_labled , subset = Patient %in% pre_NR_patient )
cell_ratio_sample_total <- prop.table(table( BCC_onlyR $Celltype..original. ,BCC_onlyR $TimePoint), margin = 2)

dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                         cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_NR_total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


# ## 2) C1QC / FCN1 TAM ratio
# tam_only_metadata <- BCC_sc_TAM_labled@meta.data[BCC_sc_TAM_labled@meta.data$Celltype..original. %in% c("FCN1+ TAM","C1QC+ TAM"),]
# cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..original.,tam_only_metadata$Treatment), margin = 2)
# 
# dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
#                              "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
#                              cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
# pdf("BCC_NR_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
# ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
#   geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
# dev.off()


## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8_ex_T_cells",],
                               cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("BCC_NR_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()

######### 14. R pre post MIF CD74 expression vs 信号强度 （archive) #####
### only CD8 effector
pre_BCC_sc_TAM_labled$Response <- ifelse(pre_BCC_sc_TAM_labled$Patient %in% pre_NR_patient,"NR", "R")
cd8_tex_only <- subset(pre_BCC_sc_TAM_labled, subset = Celltype..original. == "CD8_ex_T_cells")
expr_cd8 <- AverageExpression(cd8_tex_only, assays = "RNA", slot = "data",group.by = "Response")
expr_cd8 <- expr_cd8$RNA["MIF",]

### only C1QC TAM
c1qc_only <- subset(pre_BCC_sc_TAM_labled, subset = Celltype..original. == "C1QC+ TAM")
expr_c1qc_only <- AverageExpression(c1qc_only, assays = "RNA", slot = "data",group.by = "Response")
expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]
expr_c1qc_only <- c(expr_c1qc_only,0)
### plot
dotplot_mif_cd74 <- data.frame("MIF_expression" = expr_cd8, "CD74_expression" = expr_c1qc_only,
                               cd74_prob = c(pre_NR_mif_cd74_scores, pre_R_mif_cd74_scores),sample = names(expr_cd8))
pdf("BCC_mif_cd74_expression_vs_signal.pdf")
ggplot(dotplot_mif_cd74, aes(x=MIF_expression, y=CD74_expression, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
