library(ggpubr)
library(dplyr)
library(Seurat)

library(ggsci)
library(CellChat)
library(GSVA)
library("ggstatsplot")
library(viridis)
######### 1. Data import (In server) ########
#UVM_sc <- Read10X_h5("~/onedrive/Work/phD/phd_project/SiT/rawData/validation/Immunotherapy_metastasis /UVM_GSE141526_expression.h5")  
UVM_sc <- Read10X_h5("/public/home/ac8q5y82z9/data/0.raw/SiT_validation/validation/no_treatment/UVM_GSE139829_expression.h5")  
UVM_sc <- CreateSeuratObject(UVM_sc,project = "UVM") #后面就可以单细胞处理的标准流程啦
UVM_sc_metadata <- read.table("/public/home/ac8q5y82z9/data/0.raw/SiT_validation/validation/no_treatment/UVM_GSE139829_CellMetainfo_table.tsv",sep = "\t",row.names = 1,header = T)

######### 2. Add metadata (In server) #######
UVM_sc <- AddMetaData(
  object = UVM_sc,
  metadata = UVM_sc_metadata
)

######### 3. Scale and PCA （In server) ########
UVM_sc_feature <- FindVariableFeatures(UVM_sc)
UVM_sc_scaled <- ScaleData(UVM_sc_feature)
UVM_sc_PCA <- RunPCA(UVM_sc_scaled,features = VariableFeatures(object = UVM_sc_scaled))

UVM_sc_PCA_umap <- UVM_sc_PCA %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors( dims = 1:10) %>%
  FindClusters( dims = 1:10,resolution = 0.5) %>% 
  identity()
UVM_sc_PCA_tsne <- UVM_sc_PCA_umap %>% 
  RunTSNE( dims = 1:10)

#### save
saveRDS(UVM_sc_PCA_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_sc_PCA_tsne")

UVM_sc_PCA_tsne <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_aPD1_sc_PCA_tsne")


######### 4. Extract macrophage and cluster again （In sever) ######
mcells_only_sce_anno_ori <- subset(UVM_sc_PCA_tsne , subset = Celltype..major.lineage. == "Mono/Macro")

mcells_only_sce_anno <- SCTransform(mcells_only_sce_anno_ori) %>% RunPCA() 
pc.num=1:10
mcells_only_sce_anno <- RunUMAP(mcells_only_sce_anno, dims=pc.num) %>%
  FindNeighbors( dims=pc.num) %>% 
  FindClusters(resolution=0.8)   
mcells_only_sce_anno@active.assay = "SCT"

saveRDS(mcells_only_sce_anno,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_mcells_only_sce_anno.rds")
mcells_only_sce_anno <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_mcells_only_sce_anno.rds")

### add C1QC label
mcells_only_sce_anno$TAM_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(2,3,6,7,10,11), "FCN1+ TAM",
                                         ifelse(mcells_only_sce_anno$seurat_clusters %in% c(8), "C1QC+ TAM", "Other TAM"))

mcells_only_sce_anno$TAM_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(9), "FCN1+ TAM",
                                         ifelse(mcells_only_sce_anno$seurat_clusters %in% c(15,13), "Other TAM", "C1QC+ TAM"))

saveRDS(mcells_only_sce_anno,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_mcells_only_sce_anno_labeled.rds")
mcells_only_sce_anno <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_mcells_only_sce_anno_labeled.rds")

######### 5. Add label into original sc #######

TAM_label <- mcells_only_sce_anno@meta.data[,"TAM_label"]
names(TAM_label ) <-  rownames(mcells_only_sce_anno@meta.data)
UVM_sc_TAM_labled <- AddMetaData(
  object = UVM_sc_PCA_tsne,
  metadata = as.data.frame(TAM_label)
)

head(UVM_sc_TAM_labled@meta.data)
UVM_sc_TAM_labled$Celltype..major.lineage. <- ifelse(is.na(UVM_sc_TAM_labled$TAM_label), UVM_sc_TAM_labled$Celltype..major.lineage., UVM_sc_TAM_labled$TAM_label)
table(UVM_sc_TAM_labled$Celltype..major.lineage.)

saveRDS(UVM_sc_TAM_labled,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/new_UVM_sc_TAM_labled.rds" )

saveRDS(UVM_sc_TAM_labled,"~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled.rds" )
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/new_UVM_sc_TAM_labled.rds")


######## 6. Extract CD8 tex and cluster again （In sever) ######
cd8tex_only_sce_anno_ori <- subset(UVM_sc_TAM_labled , subset = Celltype..major.lineage. == "CD8Tex")

cd8tex_only_sce_anno <- SCTransform(cd8tex_only_sce_anno_ori) %>% RunPCA() 
pc.num=1:10
cd8tex_only_sce_anno <- RunUMAP(cd8tex_only_sce_anno, dims=pc.num) %>%
  FindNeighbors( dims=pc.num) %>% 
  FindClusters(resolution=0.8)   

cd8tex_only_sce_anno@active.assay = "SCT"

### add CD8 CXCL13 label
cd8tex_only_sce_anno$cd8tex_label <- ifelse(cd8tex_only_sce_anno$seurat_clusters %in% c(4,9), "CD8+ CXCL13+ Tex", "Other CD8Tex")

######### 7. Add label into original sc #######
cd8tex_label <- cd8tex_only_sce_anno@meta.data[,"cd8tex_label"]
names(cd8tex_label ) <-  rownames(cd8tex_only_sce_anno@meta.data)
UVM_sc_TAM_labled <- AddMetaData(
  object = UVM_sc_TAM_labled,
  metadata = as.data.frame(cd8tex_label)
)

UVM_sc_TAM_labled$Celltype..major.lineage. <- ifelse(is.na(UVM_sc_TAM_labled$cd8tex_label), UVM_sc_TAM_labled$Celltype..major.lineage., UVM_sc_TAM_labled$cd8tex_label)

table(UVM_sc_TAM_labled$Celltype..major.lineage.)
saveRDS(UVM_sc_TAM_labled,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/new_UVM_sc_TAM_labled.rds" )
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/new_UVM_sc_TAM_labled.rds")
saveRDS(UVM_sc_TAM_labled,"~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/new_UVM_sc_TAM_labled.rds" )



######### 6.1 Marker （Archive) #########
## test c1qc and fcn1 diff markers
c1qc_TAM_markers <- FindMarkers(UVM_sc_TAM_labled, group.by = "Celltype..original.",ident.1 = "C1QC+ TAM", ident.2 = "FCN1+ TAM")
c1qc_TAM_markers$gene <- rownames(c1qc_TAM_markers)
c1qc_TAM_markers = c1qc_TAM_markers%>% select(gene, everything()) %>% subset(p_val<0.05)
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
top50 = top_n(c1qc_TAM_markers ,n = 50, wt = abs(avg_log2FC))
write.csv(c1qc_TAM_markers, "~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/c1qc_TAM_markers.csv", row.names = F)
write.csv(top50, "~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/top50_c1qc_TAM_markers.csv", row.names = F)

## test c1qc diff other all cluster markers

######### 6.2 Marker visulization #########

## 1）UVM_sc_PCA_tsne
pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM/UVM_marker_plot.pdf",width = 10)
DimPlot(UVM_sc_PCA_tsne, reduction = "tsne", group.by = "Celltype..major.lineage.")
FeaturePlot(UVM_sc_PCA_tsne, reduction = "tsne",features = c("CXCL13", "C1QC", "FCN1"))
VlnPlot(UVM_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Celltype..major.lineage.")
VlnPlot(UVM_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),pt.size = 0,group.by ="Celltype..minor.lineage.")
DotPlot(UVM_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),group.by = "Celltype..major.lineage.")
DotPlot(UVM_sc_PCA_tsne,features =  c("CXCL13", "C1QC", "FCN1"),group.by = "Stage")
dev.off()

### 2）mcells_only_sce_anno
pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/UVM_mcell_marker_plot.pdf",width = 8)
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "seurat_clusters", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Celltype..major.lineage.", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "TAM_label", label = F) 
DimPlot(mcells_only_sce_anno, reduction = "tsne", group.by = "Stage", label = F) 

FeaturePlot(mcells_only_sce_anno, reduction = "tsne",features = c("C1QC", "FCN1", "SIRPA","LYZ","CD74"))
VlnPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0,group.by ="seurat_clusters")
VlnPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0,group.by ="Celltype..major.lineage.")
VlnPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0,group.by ="TAM_label")
VlnPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0,group.by ="Stage")
DotPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),group.by = "seurat_clusters")
DotPlot(mcells_only_sce_anno,features =  c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),group.by = "TAM_label")

dev.off()

DotPlot(UVM_sc_TAM_labled,features = c("FCN1", "C1QC", "CD74","SIRPA"), group.by = "Celltype..major.lineage.") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_marker_dotplot.pdf",width = 8)
DimPlot(UVM_sc_TAM_labled, reduction = "tsne", group.by = "Celltype..major.lineage.", label = F) 

DotPlot(UVM_sc_TAM_labled,features = c("FCN1", "C1QC", "CXCL13"), group.by = "Celltype..major.lineage.") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
FeaturePlot(UVM_sc_TAM_labled,reduction = "tsne",features = c("FCN1", "C1QC", "CD74","CXCL13"))

dev.off()

VlnPlot(UVM_sc_TAM_labled,features = c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0, group.by = "Celltype..major.lineage.") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))

VlnPlot(UVM_sc_TAM_labled,features = c("C1QC", "FCN1", "SIRPA","LYZ","CD74"),pt.size = 0, group.by = "seurat_clusters") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))

VlnPlot(cd8tex_only_sce_anno,features = c("CD8A","CXCL13"),pt.size = 0, group.by = "cd8tex_label") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))


######### 7. Ratio ######
### (1) All celltype vs patient ####
cell_ratio_sample <- prop.table(table( UVM_sc_TAM_labled$Celltype..major.lineage., UVM_sc_TAM_labled$Patient), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled_sample_all_ratio.pdf",width = 5, height = 5)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') + scale_fill_viridis(option="b", discrete = T) +
  coord_flip() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), legend.position = "bottom")  
dev.off()

### (2) All celltype vs group ####
cell_ratio_sample <- prop.table(table( UVM_sc_TAM_labled$Stage ,UVM_sc_TAM_labled$Celltype..major.lineage.), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)
#cell_ratio_sample <- cell_ratio_sample[cell_ratio_sample$Var1 %in% c("C1QC+ TAM", "CD8Tex", 'FCN1+ TAM'),]
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled_celltype_stage.pdf",width = 5, height = 4)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+ scale_fill_viridis(option="D", discrete = T,direction = -1) +
  coord_flip() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "bottom")  
dev.off()

### (3) Celltype piechart ######
library(lessR)
UVM_sc_TAM_labled_df <- as.data.frame(UVM_sc_TAM_labled@meta.data)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/celltype_ratio_piechart.pdf")
PieChart(Celltype..major.lineage., data = UVM_sc_TAM_labled_df,
         main = NULL, fill = "blues") 
dev.off()


### (1) boxplot all UVM （Stage) ####
group_stage_compare <- function(metadata, celltype){
  boxplot_df_all <- prop.table(table(metadata$Celltype..major.lineage.,metadata$Patient ), margin = 2)
  boxplot_df_all <- boxplot_df_all[celltype,]
  boxplot_df_all <- data.frame("celltype" = boxplot_df_all, "Patient" = names(boxplot_df_all))
  boxplot_df_all <- merge(boxplot_df_all, metadata@meta.data[,c("Patient","Stage")], by = "Patient") 
  boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
  
  ggboxplot(boxplot_df_all ,x = "Stage", y = "celltype" ,fill = "Stage") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
    stat_compare_means(method = "t.test") +
    ylab(celltype)
}

group_stage_compare <- function(metadata, celltype){
  boxplot_df_all <- prop.table(table(metadata$Patient,metadata$Celltype..major.lineage. ), margin = 2)
  boxplot_df_all <- boxplot_df_all[,celltype]
  boxplot_df_all <- data.frame("celltype" = boxplot_df_all, "Patient" = names(boxplot_df_all))
  boxplot_df_all <- merge(boxplot_df_all, metadata@meta.data[,c("Patient","Stage")], by = "Patient") 
  boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
  
  ggboxplot(boxplot_df_all ,x = "Stage", y = "celltype" ,fill = "Stage") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
    stat_compare_means(method = ) +
    ylab(celltype)
}

### ratio in all cells
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/ttest_response_comparison_ratio_boxplot.pdf",width = 5, height = 5)

group_stage_compare(UVM_sc_TAM_labled, "C1QC+ TAM")
group_stage_compare(UVM_sc_TAM_labled, "CD8+ CXCL13+ Tex")
group_stage_compare(UVM_sc_TAM_labled, "Other CD8Tex")
group_stage_compare(UVM_sc_TAM_labled, "FCN1+ TAM")
group_stage_compare(UVM_sc_TAM_labled, "CD8T")
group_stage_compare(UVM_sc_TAM_labled, "CD4Tconv")
group_stage_compare(UVM_sc_TAM_labled, "Other TAM")
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/new_Celltype_malignancy_response_comparison_ratio_boxplot.pdf",width = 5, height = 5)

group_stage_compare(UVM_sc_TAM_labled, "Immune cells" )
group_stage_compare(UVM_sc_TAM_labled, "Malignant cells")

dev.off()

### ratio in imm cells
UVM_sc_TAM_labled_imm <- subset(UVM_sc_TAM_labled, subset = Celltype..malignancy. == "Immune cells" )

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/new_response_comparison_ratio_imm_cells_boxplot.pdf",width = 5, height = 5)
group_stage_compare(UVM_sc_TAM_labled_imm, "C1QC+ TAM")
group_stage_compare(UVM_sc_TAM_labled_imm, "CD8Tex")
group_stage_compare(UVM_sc_TAM_labled_imm, "FCN1+ TAM")
group_stage_compare(UVM_sc_TAM_labled_imm, "CD8T")
group_stage_compare(UVM_sc_TAM_labled_imm, "CD4Tconv")
group_stage_compare(UVM_sc_TAM_labled_imm, "Other TAM")
dev.off()

### ratio in TAM
UVM_sc_TAM_labled_TAM <- subset(UVM_sc_TAM_labled, subset = Celltype..major.lineage. %in% c("C1QC+ TAM","FCN1+ TAM", "Other TAM") )

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/response_comparison_ratio_TAM_boxplot.pdf",width = 5, height = 5)
group_stage_compare(UVM_sc_TAM_labled_TAM, "C1QC+ TAM")
group_stage_compare(UVM_sc_TAM_labled_TAM, "FCN1+ TAM")
group_stage_compare(UVM_sc_TAM_labled_TAM, "Other TAM")
dev.off()

### (2) ggbarstats #####
###（3）CD74 exp vs CD8 effector ratio vs C1QC TAM ratio ######
cell_ratio_sample_total <- prop.table(table(UVM_sc_TAM_labled$Celltype..major.lineage.,UVM_sc_TAM_labled$Patient), margin = 2)

### only C1QC TAM
c1qc_only <- subset(UVM_sc_TAM_labled, subset = Celltype..major.lineage. == "C1QC+ TAM")
expr_c1qc_only <- AverageExpression(c1qc_only, assays = "RNA", slot = "data",group.by = "Patient")
expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]
#expr_c1qc_only <- expr_c1qc_only[]

dotplot_5 <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                        "CD8_Effctor_T" = cell_ratio_sample_total["CD8+ CXCL13+ Tex",],
                        cd74_expression = expr_c1qc_only,sample = colnames(cell_ratio_sample_total))
## add stage
stage <- UVM_sc_TAM_labled@meta.data[,c("Patient","Stage")]
stage <- stage[!duplicated(stage),]
dotplot_5 <- merge(dotplot_5, stage, by.x = "sample", by.y = "Patient")

p5 <- ggplot(dotplot_5, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_expression, color=Stage)) +
  geom_point() + scale_size(range = c(1, 10)) + scale_color_viridis(discrete = T) + theme_classic() + theme(legend.position = "top")


pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_ratio_cd74_expression_point_plot.pdf",width= 8) 
p5
dev.off()

### Archive #####

##  only TAM
UVM_sc_TAM_labled_TAM <-  UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled$Celltype..original. %in% c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_sc_TAM_labled_TAM_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(UVM_sc_TAM_labled_TAM, Response, Celltype..original., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(UVM_sc_TAM_labled_TAM,  Celltype..original., Response,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()

## CD8+ T cells
UVM_sc_TAM_labled_CD8 <-  UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled$Celltype..original. %in% c("CD8_act_T_cells", "CD8_ex_T_cells", 'CD8_mem_T_cells'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_sc_TAM_labled_CD8_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(UVM_sc_TAM_labled_CD8, Response, Celltype..original., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(UVM_sc_TAM_labled_CD8,  Celltype..original., Response,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()


### (6) boxplot Primary vs Metastasis in TAM
UVM_sc_TAM_labled_metastasis <- UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled@meta.data$ Celltype..major.lineage. %in% c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM'),]
boxplot_df_metastasis <- prop.table(table(UVM_sc_TAM_labled_metastasis$ Celltype..major.lineage., UVM_sc_TAM_labled_metastasis$Patient), margin = 2)
boxplot_df_metastasis <- boxplot_df_metastasis["C1QC+ TAM",]
boxplot_df_metastasis <- data.frame("C1QC_TAM_ratio" = boxplot_df_metastasis, "Patient" = names(boxplot_df_metastasis))

boxplot_df_metastasis <- merge(boxplot_df_metastasis, UVM_sc_TAM_labled@meta.data[,c("Patient", "Stage", "Gender")], by = "Patient" )
boxplot_df_metastasis <- boxplot_df_metastasis[!duplicated(boxplot_df_metastasis),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_C1QC_TAM_ratio_only_TAM_boxplot.pdf",width = 5, height = 5)

ggboxplot(boxplot_df_metastasis ,x = "Stage", y = "C1QC_TAM_ratio",fill = 'Stage' ) + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 


dev.off()

### (7) boxplot Primary vs Metastasis in ALL
#UVM_sc_TAM_labled_metastasis <- UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled@meta.data$ Celltype..major.lineage. %in% c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM'),]
boxplot_df_metastasis <- prop.table(table(UVM_sc_TAM_labled$ Celltype..major.lineage., UVM_sc_TAM_labled$Patient), margin = 2)
boxplot_df_metastasis <- boxplot_df_metastasis[c("C1QC+ TAM",  "CD8Tex", 'FCN1+ TAM'),]
boxplot_df_metastasis <- as.data.frame(t(boxplot_df_metastasis))
colnames(boxplot_df_metastasis) <- c("Patient", "Celltype", "Freq")

boxplot_df_metastasis <- merge(boxplot_df_metastasis, UVM_sc_TAM_labled@meta.data[,c("Patient", "Stage", "Gender")], by.x = "Patient", by.y = "Patient" )
boxplot_df_metastasis <- boxplot_df_metastasis[!duplicated(boxplot_df_metastasis),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_ratio_all_celltype_boxplot.pdf",width = 5, height = 5)

ggboxplot(boxplot_df_metastasis[boxplot_df_metastasis$Celltype == "C1QC+ TAM",] ,x = "Stage", y = "Freq",fill = 'Stage' ) + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")  + ylab("C1QC+ TAM ratio")
ggboxplot(boxplot_df_metastasis[boxplot_df_metastasis$Celltype == "FCN1+ TAM",] ,x = "Stage", y = "Freq",fill = 'Stage' ) + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") + ylab("FCN1+ TAM ratio")
ggboxplot(boxplot_df_metastasis[boxplot_df_metastasis$Celltype == "CD8Tex",] ,x = "Stage", y = "Freq",fill = 'Stage' ) + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") + ylab("CD8 Tex ratio")

dev.off()


### (8) UVM ggbarstats 

##  only TAM
UVM_sc_TAM_labled_TAM <-  UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled$Celltype..major.lineage. %in% c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_sc_TAM_labled_TAM_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(UVM_sc_TAM_labled_TAM, Stage, Celltype..major.lineage., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(UVM_sc_TAM_labled_TAM,  Celltype..major.lineage., Stage,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()

## CD8+ T cells
UVM_sc_TAM_labled_CD8 <-  UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled$Celltype..major.lineage. %in% c("CD4Tconv", "CD8T", 'CD8Tex'),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_sc_TAM_labled_CD8_ratio_ggbarstats.pdf",width = 5, height = 5)

ggbarstats(UVM_sc_TAM_labled_CD8, Stage, Celltype..major.lineage., palette = 'Set2',ggstatsplot.layer = FALSE)
ggbarstats(UVM_sc_TAM_labled_CD8,  Celltype..major.lineage., Stage,palette = 'Set3',ggstatsplot.layer = FALSE)
dev.off()


### (4) CD74 expression boxplot ######
#expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]


boxplot_6 <- data.frame("CD74_expression_in_C1QC_TAM" = expr_c1qc_only,sample = names(expr_c1qc_only))
boxplot_6 <- merge(boxplot_6, stage, by.x = "sample", by.y = "Patient")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_cd74_expression_boxplot.pdf")
ggboxplot(boxplot_6 ,x = "Stage", y = "CD74_expression_in_C1QC_TAM" ,fill = "Stage") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()


########## 12. UVM All ##########
###  extract only TAM and CD8+ Tex
#TAM_Tex_sce <- subset(UVM_sc_TAM_labled , subset = Celltype..major.lineage. == c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM',  "CD8Tex"))

# compare ligand and receptors before and after

library(CellChat)
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/cellchat_all/")

## 10.1 创建cellchat对象
cellchat <- createCellChat(UVM_sc_TAM_labled , group.by = "Celltype..major.lineage.")
cellchat <- setIdent(cellchat, ident.use = "Celltype..major.lineage.") # set "labels" as default cell identity

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

saveRDS(cellchat, file = "UVM_cellchat.rds")
cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_cellchat.rds")
#cellchat <- readRDS("UVM_cellchat.rds")
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
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ CXCL13+ Tex" & c1qc_cd8_scores$target == "C1QC+ TAM",]

pdf("c1qc_cd8_interaction.pdf",width = 10)
ggdotchart(c1qc_cd8_scores, y = "prob", x ="interaction_name_2", color = "interaction_name_2",dot.size = 6,
           add = "segments") + scale_color_viridis(discrete = T) + theme(legend.position = "none",axis.text.x = element_text(angle = 35)) 

dev.off()
########## 13. UVM Sep Cellchat （Archive)########
###  extract only TAM and CD8+ Tex
#TAM_Tex_sce <- subset(UVM_sc_TAM_labled , subset = Celltype..major.lineage. == c("C1QC+ TAM", "Other TAM", 'FCN1+ TAM',  "CD8Tex"))

pri_UVM_sce <- subset(UVM_sc_TAM_labled, subset = Stage == "Primary")
meta_UVM_sce <- subset(UVM_sc_TAM_labled, subset = Stage == "Metastatic")

setwd("sep")

### 11.1 create cellchat 
pri_cellchat <- createCellChat(pri_UVM_sce, group.by = "Celltype..major.lineage.")
meta_UVM_cellchat <- createCellChat(meta_UVM_sce, group.by = "Celltype..major.lineage.")

### 11.2 构建细胞通讯网络
## before
pri_cellchat@DB <- CellChatDB.human
pri_cellchat <- subsetData(pri_cellchat) # This step is necessary even if using the whole database
pri_cellchat <- identifyOverExpressedGenes(pri_cellchat)
pri_cellchat <- identifyOverExpressedInteractions(pri_cellchat)
pri_cellchat <- projectData(pri_cellchat, PPI.human)

pri_cellchat <- computeCommunProb(pri_cellchat)
pri_cellchat <- filterCommunication(pri_cellchat, min.cells = 10)
pri_cellchat <- computeCommunProbPathway(pri_cellchat)
pri_cellchat <- aggregateNet(pri_cellchat)

saveRDS(pri_cellchat, "UVM_pri_cellchat.rds")
pri_cellchat  <- readRDS("UVM_pri_cellchat.rds")
## after
meta_UVM_cellchat@DB <- CellChatDB.human
meta_UVM_cellchat <- subsetData(meta_UVM_cellchat) # This step is necessary even if using the whole database
meta_UVM_cellchat <- identifyOverExpressedGenes(meta_UVM_cellchat)
meta_UVM_cellchat <- identifyOverExpressedInteractions(meta_UVM_cellchat)

meta_UVM_cellchat <- computeCommunProb(meta_UVM_cellchat)
meta_UVM_cellchat <- filterCommunication(meta_UVM_cellchat, min.cells = 10)
meta_UVM_cellchat <- computeCommunProbPathway(meta_UVM_cellchat)
meta_UVM_cellchat <- aggregateNet(meta_UVM_cellchat)

saveRDS(meta_UVM_cellchat, "meta_UVM_cellchat.rds")
meta_UVM_cellchat <- readRDS("meta_UVM_cellchat.rds")

### 11.3 Merge CellChat object
total_cellchat_list <- list(Primary = pri_cellchat, Metastatic = meta_UVM_cellchat)
total_merged_cellchat <- mergeCellChat(total_cellchat_list, add.names = names(total_cellchat_list),
                                       cell.prefix = T)

saveRDS(total_merged_cellchat,"total_merged_cellchat.rds")

### 11.4 可视化

# （1） 通讯数量与强度对比

gg1 <- compareInteractions(total_merged_cellchat, show.legend = F,group = c(1,2))
gg2 <- compareInteractions(total_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("total_compareInteractions.pdf")
gg1 + gg2
dev.off()
# (2) 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_groups_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(total_merged_cellchat, weight.scale = T, measure = "weight")

dev.off()

pdf("TAM_CD8_tex_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(total_merged_cellchat, weight.scale = T, measure = "weight", sources.use = c(5), targets.use = c(2,7,9))

dev.off()

# (3) Heatmap
pdf("all_netVisual_heatmap.pdf")
netVisual_heatmap(total_merged_cellchat, measure = "weight")
dev.off()

# (4) 特定信号通路的对比
weight.max <- getMaxWeight(total_cellchat_list , attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_netVisual_circle.pdf")
for (i in 1:length(total_cellchat_list )) {
  netVisual_circle(total_cellchat_list [[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", 
                                                            names(total_cellchat_list )[i]))
}
dev.off()

pdf("MIF_netVisual_circle.pdf")
for (i in 1:length(total_cellchat_list )) {
  netVisual_aggregate(total_cellchat_list[[i]], layout = "circle", signaling = c("MIF"),
                      label.edge= F, edge.weight.max = weight.max[2], 
                      edge.width.max = 12, signaling.name = paste0("MIF  - ", 
                                                                   names(total_cellchat_list )[i]))
}
dev.off()

# (5) 通路信号强度对比分析
pdf("pri_meta_rankNet.pdf",height = 10)
rankNet(total_merged_cellchat, mode= "comparison",stacked = T, comparison = c(1,2))
dev.off()
# (6) 配体-受体对比分析
pdf("pri_meta_TAM_cd8_total_netVisual_bubble.pdf")
netVisual_bubble(total_merged_cellchat,  sources.use = c(5), targets.use = c(2,7,9),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Primary", max.dataset = 1)
netVisual_bubble(total_merged_cellchat,  sources.use = c(5), targets.use = c(2,7,9),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in Metastatic", max.dataset = 2)

dev.off()

########## 14. UVM inflammation score #######
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/scores/")
inflammation_factor <- c("TNF", "IL1A", "IL6", "CXCL8", "IL12RB1")
anti_inflammation_factor <- c("IL10", "TGFB1")
chemokines <- c("CCL2", "CCL5", "CXCL12")

TAM_related_list <- list(inflammation_factor  = inflammation_factor ,
                         anti_inflammation_factor =anti_inflammation_factor,
                         chemokines =chemokines)
## c1qc TAM
UVM_sc_TAM_labled_c1qc <- subset(UVM_sc_TAM_labled,subset = Celltype..original. == "C1QC+ TAM")
expr <- AverageExpression(UVM_sc_TAM_labled_c1qc,group.by = c("Sample"), assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

gsva.res_UVM_c1qc <- gsva( expr, TAM_related_list, method = "gsva")
gsva.res_UVM_c1qc_boxplot <- as.data.frame(t(gsva.res_UVM_c1qc))
gsva.res_UVM_c1qc_boxplot$Sample <- rownames(gsva.res_UVM_c1qc_boxplot)
gsva.res_UVM_c1qc_boxplot$celltype <- "C1QC+ TAM"
## FCN1 TAM
UVM_sc_TAM_labled_fcn1 <- subset(UVM_sc_TAM_labled,subset = Celltype..original. == "FCN1+ TAM")
expr <- AverageExpression(UVM_sc_TAM_labled_fcn1,group.by = c("Sample"), assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

gsva.res_UVM_fcn1 <- gsva( expr, TAM_related_list, method = "gsva")
gsva.res_UVM_fcn1_boxplot <- as.data.frame(t(gsva.res_UVM_fcn1))
gsva.res_UVM_fcn1_boxplot$Sample <- rownames(gsva.res_UVM_fcn1_boxplot)
gsva.res_UVM_fcn1_boxplot$celltype <- "FCN1+ TAM"

boxplot_df_c1qc <- merge(gsva.res_UVM_c1qc_boxplot, UVM_sc_TAM_labled_c1qc@meta.data[,c("Sample","Response","TimePoint","Patient")],by = "Sample")
boxplot_df_c1qc <- boxplot_df_c1qc [!duplicated(boxplot_df_c1qc ),]
boxplot_df_fcn1 <- merge(gsva.res_UVM_fcn1_boxplot, UVM_sc_TAM_labled_c1qc@meta.data[,c("Sample","Response","TimePoint","Patient")],by = "Sample")
boxplot_df_fcn1 <- boxplot_df_fcn1[!duplicated(boxplot_df_fcn1),]

boxplot_df <- rbind(boxplot_df_c1qc,boxplot_df_fcn1)

### 1) total
pdf("UVM_total_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df ,x = "celltype", y = "inflammation_factor", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df ,x = "celltype", y = "anti_inflammation_factor", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df ,x = "celltype", y = "chemokines", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

dev.off()

### 2) C1qc pre R vs NR
boxplot_df_c1qc_pre <- boxplot_df_c1qc[boxplot_df_c1qc$TimePoint == "pre",]
boxplot_df_c1qc_pre$Response <- ifelse(boxplot_df_c1qc_pre$Patient %in% pre_R_patient,"R", "NR")

pdf("UVM_c1qc_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_c1qc_pre  ,x = "Response", y = "inflammation_factor", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_pre ,x = "Response", y = "anti_inflammation_factor", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_pre ,x = "Response", y = "chemokines", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()

## 3) C1qc  R pre vs post 
boxplot_df_c1qc_R <- boxplot_df_c1qc[boxplot_df_c1qc$Patient %in% pre_R_patient,]

pdf("UVM_c1qc_R_pre_post_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_c1qc_R ,x = "TimePoint", y = "inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_R ,x = "TimePoint", y = "anti_inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_R ,x = "TimePoint", y = "chemokines", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()

## 4) C1qc  NR pre vs post 
boxplot_df_c1qc_NR <- boxplot_df_c1qc[boxplot_df_c1qc$Patient %in% pre_NR_patient,]

pdf("UVM_c1qc_NR_pre_post_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_c1qc_NR ,x = "TimePoint", y = "inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_NR ,x = "TimePoint", y = "anti_inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_c1qc_NR ,x = "TimePoint", y = "chemokines", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()
### 5) FCN1 pre R vs NR
boxplot_df_fcn1_pre <- boxplot_df_fcn1[boxplot_df_fcn1$TimePoint == "pre",]
boxplot_df_fcn1_pre$Response <- ifelse(boxplot_df_fcn1_pre$Patient %in% pre_R_patient,"R", "NR")

pdf("UVM_fcn1_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_fcn1_pre  ,x = "Response", y = "inflammation_factor", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_pre ,x = "Response", y = "anti_inflammation_factor", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_pre ,x = "Response", y = "chemokines", color = "Response") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()

## 6) fcn1 R pre vs post 
boxplot_df_fcn1_R <- boxplot_df_fcn1[boxplot_df_fcn1$Patient %in% pre_R_patient,]

pdf("UVM_fcn1_R_pre_post_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_fcn1_R ,x = "TimePoint", y = "inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_R ,x = "TimePoint", y = "anti_inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_R ,x = "TimePoint", y = "chemokines", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()

## 7) fcn1  NR pre vs post 
boxplot_df_fcn1_NR <- boxplot_df_fcn1[boxplot_df_fcn1$Patient %in% pre_NR_patient,]

pdf("UVM_fcn1_NR_pre_post_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_fcn1_NR ,x = "TimePoint", y = "inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_NR ,x = "TimePoint", y = "anti_inflammation_factor", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

ggboxplot(boxplot_df_fcn1_NR ,x = "TimePoint", y = "chemokines", color = "TimePoint") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means() 

dev.off()

########## 15. UVM inflammation score #######
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/scores/")
inflammation_factor <- c("TNF", "IL1A", "IL6", "CXCL8", "IL12RB1")
anti_inflammation_factor <- c("IL10", "TGFB1","IL4","CSF1")
chemokines <- c("CCL2", "CCL5", "CXCL12")

TAM_related_list <- list(inflammation_factor  = inflammation_factor ,
                         anti_inflammation_factor =anti_inflammation_factor,
                         chemokines =chemokines)
## c1qc TAM
UVM_sc_TAM_labled_c1qc <- subset(UVM_sc_TAM_labled,subset = Celltype..major.lineage. == "C1QC+ TAM")
expr <- AverageExpression(UVM_sc_TAM_labled_c1qc,group.by = c("Patient"), assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

gsva.res_UVM_c1qc <- gsva( expr, TAM_related_list, method = "gsva")
gsva.res_UVM_c1qc_boxplot <- as.data.frame(t(gsva.res_UVM_c1qc))
gsva.res_UVM_c1qc_boxplot$Patient <- rownames(gsva.res_UVM_c1qc_boxplot)
gsva.res_UVM_c1qc_boxplot$celltype <- "C1QC+ TAM"
## FCN1 TAM
UVM_sc_TAM_labled_fcn1 <- subset(UVM_sc_TAM_labled,subset = Celltype..major.lineage. == "FCN1+ TAM")
expr <- AverageExpression(UVM_sc_TAM_labled_fcn1,group.by = c("Patient"), assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

gsva.res_UVM_fcn1 <- gsva( expr, TAM_related_list, method = "gsva")
gsva.res_UVM_fcn1_boxplot <- as.data.frame(t(gsva.res_UVM_fcn1))
gsva.res_UVM_fcn1_boxplot$Patient <- rownames(gsva.res_UVM_fcn1_boxplot)
gsva.res_UVM_fcn1_boxplot$celltype <- "FCN1+ TAM"

boxplot_df_c1qc <- merge(gsva.res_UVM_c1qc_boxplot, UVM_sc_TAM_labled_c1qc@meta.data[,c("Stage","Patient")],by = "Patient")
boxplot_df_c1qc <- boxplot_df_c1qc [!duplicated(boxplot_df_c1qc ),]
boxplot_df_fcn1 <- merge(gsva.res_UVM_fcn1_boxplot, UVM_sc_TAM_labled_c1qc@meta.data[,c("Stage","Patient")],by = "Patient")
boxplot_df_fcn1 <- boxplot_df_fcn1[!duplicated(boxplot_df_fcn1),]

boxplot_df <- rbind(boxplot_df_c1qc,boxplot_df_fcn1)

### 1) total
pdf("UVM_total_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df ,x = "celltype", y = "inflammation_factor", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df ,x = "celltype", y = "anti_inflammation_factor", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df ,x = "celltype", y = "chemokines", color = "celltype") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

dev.off()

### 2) C1qc Pri vs Meta
pdf("UVM_C1QC_pre_meta_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_c1qc ,x = "Stage", y = "inflammation_factor", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df_c1qc  ,x = "Stage", y = "anti_inflammation_factor", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df_c1qc  ,x = "Stage", y = "chemokines", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

dev.off()


### 3) FCN1 Pri vs Meta
pdf("UVM_FCN1_pre_meta_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(boxplot_df_fcn1 ,x = "Stage", y = "inflammation_factor", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df_fcn1  ,x = "Stage", y = "anti_inflammation_factor", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(boxplot_df_fcn1  ,x = "Stage", y = "chemokines", color = "Stage") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

dev.off()

######## 16. UVM MIF-CD74和 TAM  CD8 CXCL13 ratio correlation ########
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/sep")
pri_cellchat  <- readRDS("UVM_pri_cellchat.rds")
meta_cellchat <- readRDS("meta_UVM_cellchat.rds")
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled.rds")

### pri
pri_mif_cd74_scores <- subsetCommunication(pri_cellchat,slot.name = "netP")
pri_mif_cd74_scores <- pri_mif_cd74_scores[pri_mif_cd74_scores$source == "CD8Tex" & pri_mif_cd74_scores$target == "C1QC+ TAM",]
pri_mif_cd74_scores <- pri_mif_cd74_scores[pri_mif_cd74_scores$pathway_name == "MIF","prob"]

### meta
meta_mif_cd74_scores <- subsetCommunication(meta_cellchat,slot.name = "netP")
meta_mif_cd74_scores <- meta_mif_cd74_scores[meta_mif_cd74_scores$source == "CD8Tex" & meta_mif_cd74_scores$target == "C1QC+ TAM",]
meta_mif_cd74_scores <- meta_mif_cd74_scores[meta_mif_cd74_scores$pathway_name == "MIF","prob"]

## 1) total ratio
cell_ratio_sample_total <- prop.table(table( UVM_sc_TAM_labled$Celltype..major.lineage. ,UVM_sc_TAM_labled$Stage), margin = 2)

dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
                         cd74_prob = c(meta_mif_cd74_scores ,pri_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


## 2) C1QC / FCN1 TAM ratio
tam_only_metadata <- UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled@meta.data$Celltype..major.lineage. %in% c("FCN1+ TAM","C1QC+ TAM"),]
cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..major.lineage.,tam_only_metadata$Stage), margin = 2)

dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
                             "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
                             cd74_prob = c(meta_mif_cd74_scores, pri_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_C1QC_TAM_FCN1_ratio_CD8_effector.pdf")
ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()

## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8Tex",],
                             cd74_prob = c(meta_mif_cd74_scores, pri_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_dotplot_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()

######## 17. UVM pre post NR MIF-CD74和 TAM  CD8 CXCL13 ratio correlation ########
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_pre_post/")
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled.rds")

post_NR_cellchat  <- readRDS("post_NR_cellchat.rds")
pre_NR_cellchat <- readRDS("pre_NR_cellchat.rds")

### post_NR
post_NR_mif_cd74_scores <- subsetCommunication(post_NR_cellchat,slot.name = "netP")
post_NR_mif_cd74_scores <- post_NR_mif_cd74_scores[post_NR_mif_cd74_scores$source == "CD8Tex" & post_NR_mif_cd74_scores$target == "C1QC+ TAM",]
post_NR_mif_cd74_scores <- post_NR_mif_cd74_scores[post_NR_mif_cd74_scores$pathway_name == "MIF","prob"]

### pre_NR
pre_NR_mif_cd74_scores <- subsetCommunication(pre_NR_cellchat,slot.name = "netP")
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$source == "CD8Tex" & pre_NR_mif_cd74_scores$target == "C1QC+ TAM",]
pre_NR_mif_cd74_scores <- pre_NR_mif_cd74_scores[pre_NR_mif_cd74_scores$pathway_name == "MIF","prob"]

## 1) total ratio
cell_ratio_sample_total <- prop.table(table( UVM_sc_TAM_labled$Celltype..original. ,UVM_sc_TAM_labled$Treatment), margin = 2)

dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
                         cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_NR_total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


## 2) C1QC / FCN1 TAM ratio
tam_only_metadata <- UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled@meta.data$Celltype..original. %in% c("FCN1+ TAM","C1QC+ TAM"),]
cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..original.,tam_only_metadata$Treatment), margin = 2)

dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
                             "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
                             cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_NR_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()


## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8Tex",],
                               cd74_prob = c(post_NR_mif_cd74_scores ,pre_NR_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_NR_dotplot_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()


######## 18. UVM pre post R MIF-CD74和 TAM  CD8 CXCL13 ratio correlation ########
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_pre_post/")
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled.rds")

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

## 1) total ratio
UVM_onlyR <- subset(UVM_sc_TAM_labled , subset = Patient %in% pre_R_patient )
cell_ratio_sample_total <- prop.table(table( UVM_onlyR $Celltype..original. ,UVM_onlyR $Treatment), margin = 2)

dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8_ex_T_cells",],
                         cd74_prob = c(post_R_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_R_total_ratio_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
#### MIF CD74 跟CD8 effector ratio  c1qc tam 都反着 有点子不好解释 但是可以说cd8 effector和 c1qc tam共定位导致的


## 2) C1QC / FCN1 TAM ratio
tam_only_metadata <- UVM_sc_TAM_labled@meta.data[UVM_sc_TAM_labled@meta.data$Celltype..original. %in% c("FCN1+ TAM","C1QC+ TAM"),]
cell_ratio_sample_TAM <- prop.table(table( tam_only_metadata$Celltype..original.,tam_only_metadata$Treatment), margin = 2)

dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
                             "CD8_Effctor_T" = cell_ratio_sample_total["CD8Tex",],
                             cd74_prob = c(post_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_NR_dotplot_C1QC_TAM_ratio_CD8_effector.pdf")
ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()


## 3）C1QC / Cd8  TAM ratio

dotplot_df_ratio <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8Tex",],
                               cd74_prob = c(post_NR_mif_cd74_scores ,pre_R_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
pdf("UVM_NR_dotplot_C1QC_TAM_CD8_effector_ratio.pdf")
ggplot(dotplot_df_ratio, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()


#### 13.  UVM MIF CD74 expression vs 信号强度 #####
UVM_sc_TAM_labled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_sc_TAM_labled.rds")

### only CD8 effector
cd8_tex_only <- subset(UVM_sc_TAM_labled, subset = Celltype..major.lineage. == "CD8Tex")
expr_cd8 <- AverageExpression(cd8_tex_only, assays = "RNA", slot = "data",group.by = "Stage")
expr_cd8 <- expr_cd8$RNA["MIF",]

### only C1QC TAM
c1qc_only <- subset(UVM_sc_TAM_labled, subset = Celltype..major.lineage. == "C1QC+ TAM")
expr_c1qc_only <- AverageExpression(c1qc_only, assays = "RNA", slot = "data",group.by = "Stage")
expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]

### plot
dotplot_mif_cd74 <- data.frame("MIF_expression" = expr_cd8, "CD74_expression" = expr_c1qc_only,
                               cd74_prob = c(meta_mif_cd74_scores ,pri_mif_cd74_scores),sample = names(expr_cd8))
pdf("UVM_mif_cd74_expression_vs_signal.pdf")
ggplot(dotplot_mif_cd74, aes(x=MIF_expression, y=CD74_expression, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(5, 15)) + scale_color_viridis(discrete = T)
dev.off()
