######### 0. Library import ########
library(Seurat)
library(GSVA)
library(survival)
library(ggfortify)
library(survminer)
library(readxl)
library(ggstatsplot)
library(harmony)
library(viridis)
library(tibble)
library(VennDiagram)
######## 1.Read data and creat SeuratObject #####
# GSE203067
LN_validataion <- CreateSeuratObject(counts = Read10X("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/"), 
                   project = "GSE203067" )

LN_validataion@meta.data$sample <- sapply(strsplit(rownames(LN_validataion@meta.data),"-"),"[",2)
table(LN_validataion@meta.data$group )
## keep all
#LN_validataion_ <- subset(LN_validataion, group %in% c(1,2,3,4))
#table(LN_validataion_@meta.data$group )
#LN_validataion_ori <- LN_validataion
table(LN_validataion$group)
saveRDS(LN_validataion,"/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067//LN_validataion.rds")
LN_validataion <- readRDS("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067//LN_validataion.rds")
LN_validataion <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/LN_validataion.rds")

## addmeta
metadata <- read.delim("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/ESCC_Jia.metadata.txt")
rownames(metadata) <- metadata$cells
nrow(metadata)
table(metadata$group)
LN_validataion$cell_names <- rownames(LN_validataion@meta.data)
LN_validataion <- subset(LN_validataion, subset = cell_names %in% metadata$cells )
LN_validataion<- AddMetaData(LN_validataion, metadata =metadata)
head(LN_validataion@meta.data)
table(LN_validataion$sample)

### 2. Rclustering results ####
# LN_validataion_qc <- NormalizeData(LN_validataion_anno, normalization.method = "LogNormalize") 
# LN_validataion_qc <- FindVariableFeatures(LN_validataion_qc,nfeatures = 2000)
# all.genes <- rownames(LN_validataion_qc )
# LN_validataion_qc <- ScaleData(LN_validataion_qc,features = all.genes)
# LN_validataion_normalized <- RunPCA(LN_validataion_qc, features = VariableFeatures(object = LN_validataion_qc))
# 
# LN_validataion <- RunHarmony(LN_validataion_normalized, group.by.vars = "group_sample") 
# LN_validataion_harmony_umap <- LN_validataion_harmony %>% 
#   RunUMAP(reduction = "harmony", dims = 1:10,n.neighbors = 50,n.epochs=200,negative.sample.rate = 50) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
#   FindClusters( dims = 1:0,resolution = seq(0.2,0.8,0.2)) %>% 
#   identity()

LN_validataion <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/LN_validataion.rds")
LN_validataion_SCT <- SCTransform(object = LN_validataion, assay = "RNA",ncells=5000,conserve.memory=TRUE)
LN_validataion_PCA <- RunPCA(LN_validataion_SCT, features = VariableFeatures(object = LN_validataion_SCT))
LN_validataion_umap <- LN_validataion_PCA %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors( dims = 1:30) %>%
  FindClusters( dims = 1:30,resolution = 0.1) %>%
  identity()
LN_validataion_tsne <- LN_validataion_umap %>%
  RunTSNE( dims = 1:10)
saveRDS(LN_validataion_SCT,"/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/LN_validataion_SCT.rds")
saveRDS(LN_validataion_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/LN_validataion_tsne.rds")


LN_validataion_tsne <- readRDS("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/LN_validataion_tsne.rds")

##### 3. Dimplot ######
#pdf("/public/home/ac8q5y82z9/data/3.Results/LY_scRNA/fibs/fib_dimplot_pc20_0.5.pdf")
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/GSE/dimplot_pc10_0_1.pdf")
DimPlot(LN_validataion_tsne, reduction = "tsne", group.by = "cell_type", pt.size=0.8, raster = F, label = T)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_tsne, reduction = "tsne", group.by = "sample", pt.size=0.5, raster = F, label = T)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_tsne, reduction = "umap", group.by = "cell_type", pt.size=0.5, raster = F, label = T)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_tsne, reduction = "umap", group.by = "sample", pt.size=0.5, raster = F, label = T)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/GSE_ESCC_mportant_genes_Featureplot.pdf")
#pdf()
FeaturePlot(LN_validataion_tsne, features = c("FCN1","CD74","C1QC","CXCL13"))
dev.off()
####### 4. Macrophages ######
macro_LN_validation <- subset(LN_validataion_tsne, subset = cell_type == "Myeloid cells")
table(macro_LN_validation$CellSubtype)
VlnPlot(macro_LN_validation, features=c("FCN1","C1QC","CD74"), group.by = "subtypes")

macro_LN_validation$new_groups <- ifelse(macro_LN_validation$subtypes %in% c("C2_Monocytes", "C0_Monocytes"),"FCN1+ TAMs",
                                         ifelse(macro_LN_validation$subtypes %in% c("C10_Macrophages", "C1_Macrophages",
                                                                                    "C11_Macrophages","C12_Macrophages",
                                                                                    "C7_Macrophages","C4_Macrophages",
                                                                                    "C3_Macrophages","C6_Macrophages","C17_Macrophages"),"C1QC+ TAMs", "Other myeloids and TAMs"))
VlnPlot(macro_LN_validation, features=c("FCN1","C1QC","CD74"), group.by = "new_groups")


####### 5. T cells ####
tcell_LN_validation <- subset(LN_validataion_tsne, subset = cell_type == "T cells")
table(tcell_LN_validation$group)
table(tcell_LN_validation$CellSubtype)
table(tcell_LN_validation$subtypes)

VlnPlot(tcell_LN_validation, features=c("CXCL13"), group.by = "subtypes")
tcell_LN_validation$new_groups <- ifelse(tcell_LN_validation$subtypes %in% c("C0_EX", 
                                            "C16_EX","C21_EX","C27_FH2","C26_EX",
                                            "C18_Proliferating","C13_FH1",""),"CD8+ CXCL13+ Tcells",
                                         tcell_LN_validation$CellSubtype )
 
tcell_LN_validation$new_groups <- ifelse(tcell_LN_validation$CellSubtype %in% c("EX","FH1","FH2","Proliferating"),"CD8+ CXCL13+ Tcells",
                                         tcell_LN_validation$CellSubtype )

VlnPlot(tcell_LN_validation, features=c("CXCL13"), group.by = "new_groups")


####### 6. Re annotation data #####
LN_validataion_tsne <- LN_validataion_SCT
# LN_validataion_tsne$new_groups <- ifelse(LN_validataion_tsne$cell_type == "T cells", "other T cells",
#                                           ifelse(LN_validataion_tsne$CellSubtype %in% c("EX","FH1","FH2","Proliferating") & LN_validataion_tsne$cell_type == "T cells","CD8+ CXCL13+ Tcells",
#                                    ifelse(LN_validataion_tsne$subtypes %in% c("C2_Monocytes", "C0_Monocytes"),"FCN1+ TAMs",
#                                           ifelse(LN_validataion_tsne$subtypes %in% c("C10_Macrophages", "C1_Macrophages",
#                                                                                      "C11_Macrophages","C12_Macrophages",
#                                                                                      "C7_Macrophages","C4_Macrophages",
#                                                                                      "C3_Macrophages","C6_Macrophages","C17_Macrophages"),"C1QC+ TAMs",
#                                                  ifelse(LN_validataion_tsne$CellSubtype == "Myloid cells", "Other myeloids and TAMs", LN_validataion_tsne$cell_type)))))

LN_validataion_tsne$new_groups <- ifelse(LN_validataion_tsne$CellSubtype %in% c("EX","FH1","FH2","Proliferating") ,"CD8+ CXCL13+ Tcells",
                                                ifelse(LN_validataion_tsne$cell_type == "T cells", "other T cells",
                                                       ifelse(LN_validataion_tsne$subtypes %in% c("C2_Monocytes", "C0_Monocytes"),"FCN1+ TAMs",
                                                              ifelse(LN_validataion_tsne$subtypes %in% c("C10_Macrophages", "C1_Macrophages",
                                                                                                         "C11_Macrophages","C12_Macrophages",
                                                                                                         "C7_Macrophages","C4_Macrophages",
                                                                                                         "C3_Macrophages","C6_Macrophages","C17_Macrophages"),"C1QC+ TAMs",
                                                                     ifelse(LN_validataion_tsne$cell_type == "Myeloid cells", "Other myeloids and TAMs",LN_validataion_tsne$cell_type)))))
                                                
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/GSE/re_anno_dimplot_pc10_0_1.pdf")
pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/validation/re_anno_dimplot_pc10_0_1.pdf")
DimPlot(LN_validataion_re_anno, reduction = "umap", group.by = "new_groups", pt.size=0.8, raster = F, label = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_re_anno, reduction = "tsne", group.by = "new_groups", pt.size=0.8, raster = F, label = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_re_anno, reduction = "umap", group.by = "cell_type", pt.size=0.8, raster = F, label = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(LN_validataion_re_anno, reduction = "tsne", group.by = "cell_type", pt.size=0.8, raster = F, label = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
dev.off()

saveRDS(LN_validataion_tsne,"/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/LN_validataion_re_anno.rds")
saveRDS(LN_validataion_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/validation/LN_validataion_re_anno.rds")
LN_validataion_re_anno <- readRDS("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/GSE203067/LN_validataion_re_anno.rds")

######### 7. Ratio ######
### (1) All celltype vs patient ####
LN_validataion_tsne$sample <- sapply(strsplit(rownames(LN_validataion_tsne@meta.data),"-"),"[",2)
cell_ratio_sample <- prop.table(table( LN_validataion_tsne$new_groups, LN_validataion_tsne$sample), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_validataion_tsne_sample_all_ratio.pdf",width = 5, height = 5)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') + scale_fill_viridis(option="b", discrete = T) +
  coord_flip() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), legend.position = "bottom")  
dev.off()

### (2) All celltype vs group ####
cell_ratio_sample <- prop.table(table( LN_validataion_tsne$group ,LN_validataion_tsne$new_groups), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)
#cell_ratio_sample <- cell_ratio_sample[cell_ratio_sample$Var1 %in% c("C1QC+ TAM", "CD8Tex", 'FCN1+ TAM'),]
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_validataion_tsne_celltype_group.pdf",width = 5, height = 4)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+ scale_fill_viridis(option="A", discrete = T,direction = -1) +
  coord_flip() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "bottom")  
dev.off()

### (3) Celltype piechart ######
library(lessR)
LN_validataion_tsne_df <- as.data.frame(LN_validataion_tsne@meta.data)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_celltype_ratio_piechart.pdf")
PieChart(new_groups, data = LN_validataion_tsne_df,
         main = NULL, fill = "blues") 
dev.off()

####### 8. Dotplot #####
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_marker_plot.pdf",width = 10)
DotPlot(LN_validataion_tsne, features = c("FCN1","C1QC","CXCL13"), group.by = "new_groups") + scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
dev.off()

##### 9. Ratio boxplot ####
group_stage_compare <- function(metadata, celltype){
  boxplot_df_all <- prop.table(table(metadata$new_groups,metadata$sample ), margin = 2)
  boxplot_df_all <- boxplot_df_all[celltype,]
  boxplot_df_all <- data.frame("new_groups" = boxplot_df_all, "sample" = names(boxplot_df_all))
  boxplot_df_all <- merge(boxplot_df_all, metadata@meta.data[,c("groups","sample")], by = "sample") 
  boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
  
  ggboxplot(boxplot_df_all ,x = "group", y = "new_groups" ,fill = "group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
    stat_compare_means(method = "t.test") +
    ylab(celltype)
}

group_stage_compare <- function(metadata, celltype){
  boxplot_df_all <- prop.table(table(metadata$Patient,metadata$new_groups ), margin = 2)
  boxplot_df_all <- boxplot_df_all[,celltype]
  boxplot_df_all <- data.frame("celltype" = boxplot_df_all, "Patient" = names(boxplot_df_all))
  boxplot_df_all <- merge(boxplot_df_all, metadata@meta.data[,c("Patient","Stage")], by = "Patient") 
  boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
  
  ggboxplot(boxplot_df_all ,x = "Stage", y = "celltype" ,fill = "Stage") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
    stat_compare_means(method = ) +
    ylab(celltype)
}

### ratio in all cells
metadata <- LN_validataion_tsne@meta.data
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_group_response_comparison_ratio_boxplot.pdf",width = 5, height = 5)
boxplot_df_all <- prop.table(table(metadata$new_groups,metadata$sample ), margin = 2)
boxplot_df_all <- boxplot_df_all["C1QC+ TAMs",]
boxplot_df_all <- data.frame("new_groups" = boxplot_df_all, "sample" = names(boxplot_df_all))
boxplot_df_all <- merge(boxplot_df_all, metadata[,c("group","sample")], by = "sample") 
boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
boxplot_df_all$group <- factor(boxplot_df_all$group,c("ESCC_N","ESCC_PT","ESCC_LN","ESCC_LNM") )

p1 <- ggboxplot(boxplot_df_all ,x = "group", y = "new_groups" ,fill = "group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("ESCC_N","ESCC_PT"), c("ESCC_LN","ESCC_LNM"))) +
  ylab("C1QC+ TAMs")

boxplot_df_all <- prop.table(table(metadata$new_groups,metadata$sample ), margin = 2)
boxplot_df_all <- boxplot_df_all["CD8+ CXCL13+ Tcells",]
boxplot_df_all <- data.frame("new_groups" = boxplot_df_all, "sample" = names(boxplot_df_all))
boxplot_df_all <- merge(boxplot_df_all, metadata[,c("group","sample")], by = "sample") 
boxplot_df_all <- boxplot_df_all[!duplicated(boxplot_df_all),]
boxplot_df_all$group <- factor(boxplot_df_all$group,c("ESCC_N","ESCC_PT","ESCC_LN","ESCC_LNM") )

p2 <- ggboxplot(boxplot_df_all ,x = "group", y = "new_groups" ,fill = "group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("ESCC_N","ESCC_PT"), c("ESCC_LN","ESCC_LNM"))) +
  ylab("CD8+ CXCL13+ Tcells")
print(p1)
print(p2)

dev.off()

### (4) CD74 expression boxplot ######
#expr_c1qc_only <- expr_c1qc_only$RNA["CD74",]
c1qc_only <- subset(LN_validataion_tsne, subset = new_groups == "C1QC+ TAMs")
expr_c1qc_only <- AverageExpression(c1qc_only , assays = "SCT", slot = "data",group.by = "sample")
expr_c1qc_only <- expr_c1qc_only$SCT["CD74",]
#expr_c1qc_only <- expr_c1qc_only[expr_c1qc_only != Inf]
#expr_c1qc_only <- scale(as.numeric(expr_c1qc_only))

boxplot_6 <- data.frame("CD74_expression_in_C1QC_TAM" = expr_c1qc_only,sample = names(expr_c1qc_only))
boxplot_6 <- merge(boxplot_6, metadata[,c("group","sample")], by = "sample")
boxplot_6 <- boxplot_6[!duplicated(boxplot_6),]
boxplot_6$group <- factor(boxplot_6$group,c("ESCC_N","ESCC_PT","ESCC_LN","ESCC_LNM") )
pdf("~/onedrive/Work/phD/phd_project/SiT/results/validation/LN_cd74_expression_boxplot.pdf")
ggboxplot(boxplot_6 ,x = "group", y = "CD74_expression_in_C1QC_TAM" ,fill = "group") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("ESCC_N","ESCC_PT"), c("ESCC_LN","ESCC_LNM"))) 
dev.off()

########## 12.  cellchat All ##########

library(CellChat)
setwd("~/onedrive/Work/phD/phd_project/SiT/results/validation/ESCC/cellchat_only_LN/")

######### 13. cellchat Only LN #######
LN_validataion_tsne_only_LN <- subset(LN_validataion_tsne, subset = group %in% c("ESCC_LN", "ESCC_LNM"))

## 10.1 创建cellchat对象
#LN_validataion_tsne_ds<- subset(LN_validataion_tsne, downsample = 50000)
cellchat <- createCellChat(LN_validataion_tsne_only_LN , group.by = "new_groups")
cellchat <- setIdent(cellchat, ident.use = "new_groups") # set "labels" as default cell identity

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

saveRDS(cellchat, file = "only_LN_ESCC_cellchat.rds")
#cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/UVM_cellchat.rds")
#cellchat <- readRDS("UVM_cellchat.rds")
#### 10.3 可视化 
## (0) save csv
 

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
p1 <- netVisual_bubble(cellchat,  targets.use = c(2),
                       remove.isolate = F)
p2 <- netVisual_bubble(cellchat,  targets.use = c(3),
                       remove.isolate = F)
p3 <- netVisual_bubble(cellchat,  sources.use = c(2),
                       remove.isolate = F)
p4 <- netVisual_bubble(cellchat,  sources.use = c(3),
                       remove.isolate = F)



pdf("all_MIF_cd8_effector_all_TAM.pdf")
p1
p2
p3
p4
dev.off()

##### Lolipop
c1qc_cd8_scores <- subsetCommunication(cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ CXCL13+ Tcells" & c1qc_cd8_scores$target == "C1QC+ TAMs",]

pdf("c1qc_cd8_interaction.pdf",width = 10)
ggdotchart(c1qc_cd8_scores, y = "prob", x ="interaction_name_2", color = "interaction_name_2",dot.size = 6,
           add = "segments") + scale_color_viridis(discrete = T) + theme(legend.position = "none",axis.text.x = element_text(angle = 35)) 

dev.off()
saveRDS(cellchat, file = "only_LN_ESCC_cellchat.rds")


##### 14. Baidu II data validation #####

### Count
count <- read.table("/Users/lexiiiii/onedrive/Work/phD/phd_project/TME_gender/rawdata/baiduII_gene_count_.txt",header = T, row.names = 1)
nrow(count)
nrow(na.omit(count))
count <- count[,grep("T",colnames(count))]
colnames(count) <- gsub("X","",colnames(count))
colnames(count) <- gsub("T","",colnames(count))
row_list <- c()
for (rowname in rownames(count)) {
  row_list <- c(row_list , strsplit(rowname, "_")[[1]][2])
}
count$gene_name <- row_list

### Metadata
metadata_baiduII_ori <- read.csv(file = "/Users/lexiiiii/onedrive/Work/phD/phd_project/TME_gender/rawdata/baiduII_clinical_data.csv", header = TRUE)
metadata_baiduII_ori <- metadata_baiduII_ori[,c(1:30)]
metadata_BaiduII <- metadata_baiduII_ori[,c("sample.ID............BDESCC2..", "Gender",
                                            "Smoking.history",  "Grade","TNM.stage.the.Eighth.Edition.","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.")]
colnames(metadata_BaiduII) <- c("Tumor_Sample_Barcode","Gender","Smoking status","Grade","TNM stage","Metastasis_recurrence_status")
metadata_BaiduII$`TNM stage` <- ifelse(metadata_BaiduII$`TNM stage` %in% c("IIIA","IIIB"),3,
                                       ifelse(metadata_BaiduII$`TNM stage` %in% c("IIA","ⅡA","IIB"),2,
                                              ifelse(metadata_BaiduII$`TNM stage` == "IB",1,
                                                     ifelse(metadata_BaiduII$`TNM stage` == "ⅣA",4,"stop"))))

##### 3. Marker choice C1QC CD74 LYZ CD68 MS4A6A  ######
## Get CD8+ CXCL13+ Tex cells markers 
cd8_cxcl13_markers <- read.csv("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/res_0.5_cosg_cluster_marker.csv")
cd8_cxcl13_markers <- c(cd8_cxcl13_markers$X11,"CD8A","CXCL13","CD3E")
## Get C1QC+ CD74+ Macrophages markers
c1qc_TAM_markers <- read.csv("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/macrophage/res_2_mcells_cosg_cluster_marker.csv")
c1qc_TAM_markers <- c(c1qc_TAM_markers$X4,"CD74","LYZ")
## GSVA 
gsva_marker_list <- list(c1qc_TAM_markers,cd8_cxcl13_markers )
names(gsva_marker_list) <- c("C1QC+ CD74+ Macrophages","CD8+ CXCL13+ T cells")

count_co_interaction <- count[count$gene_name %in% c(c1qc_TAM_markers,cd8_cxcl13_markers ),]
## remove dplicated gene
count_co_interaction <- count_co_interaction[!duplicated(count_co_interaction$gene_name),]
rownames(count_co_interaction) <- count_co_interaction$gene_name
count_co_interaction <- count_co_interaction[,1:ncol( count_co_interaction)-1]

count_co_interaction <- as.matrix(count_co_interaction)
gsva.res <- gsva(count_co_interaction, gsva_marker_list, method = "gsva")
gsva.res <- t(gsva.res)
## High group vs low group 
### macrophage ##
macro_patients_high <- rownames(gsva.res[gsva.res[,1] >= mean(gsva.res[,1]),])
#macro_patients_high <- as.numeric(macro_patients_high )
macro_patients_low <- rownames(gsva.res[gsva.res[,1] < mean(gsva.res[,1]),])

### Tcells ##
tcells_patients_high <- rownames(gsva.res[gsva.res[,2] >= mean(gsva.res[,2]),])
#tcells_patients_high <- as.numeric(tcell_patients_high )
tcells_patients_low <- rownames(gsva.res[gsva.res[,2] < mean(gsva.res[,2]),])



##### 4. Survival plot #######
metadata_baiduII_survival <- metadata_baiduII_ori[,c("sample.ID............BDESCC2..", "overall.survival...............1..dead..0.alive.",
                                                     "overall.survival........time","Gender","Age", 
                                                     "Location","Family.History.of.ESCC","Family.History.of.other.cancers",
                                                     "recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.",
                                                     "TNM.stage.the.Eighth.Edition.","TNM.stage..the.7th.Edition.","Grade",
                                                     "Smoking.history","Drinking.history","T","N","M","TNM.status")]
colnames(metadata_baiduII_survival)<- c("Sample_ID","Survival_status","Survival_time",
                                        "Gender", "Age", "Location", "ESCC_family_history",
                                        "Other_cancer_family_history", "Recurrence_status",
                                        "TNM_stage_8","TNM_stage_7","Grade","smoking",
                                        "drinking","T","N","M","TNM.status")

## filter
metadata_baiduII_survival <- metadata_baiduII_survival[metadata_baiduII_survival$Survival_status != "loss",]
metadata_baiduII_survival$Survival_status <- as.numeric(metadata_baiduII_survival$Survival_status)
metadata_baiduII_survival$Survival_time <- as.numeric(metadata_baiduII_survival$Survival_time)

## Add meta
metadata_baiduII_survival$Macrophage <- ifelse(metadata_baiduII_survival$Sample_ID %in% macro_patients_high, "C1QC+ CD74+ Macrophages high","C1QC+ CD74+ Macrophages low")
metadata_baiduII_survival$Tcells <- ifelse(metadata_baiduII_survival$Sample_ID %in% tcells_patients_high, "CD8+ CXCL13+ T cells high","CD8+ CXCL13+ T cells low")
metadata_baiduII_survival$Total <- ifelse(metadata_baiduII_survival$Sample_ID %in% macro_patients_high & metadata_baiduII_survival$Sample_ID %in% tcells_patients_high, "C1QC+ CD74+ Macrophages high CD8+ CXCL13+ T cells high",
                                          ifelse(metadata_baiduII_survival$Sample_ID %in% macro_patients_high & metadata_baiduII_survival$Sample_ID %in% tcells_patients_low, "C1QC+ CD74+ Macrophages high CD8+ CXCL13+ T cells low",
                                                 ifelse(metadata_baiduII_survival$Sample_ID %in% macro_patients_low & metadata_baiduII_survival$Sample_ID %in% tcells_patients_low, "C1QC+ CD74+ Macrophages low CD8+ CXCL13+ T cells low",
                                                        ifelse(metadata_baiduII_survival$Sample_ID %in% macro_patients_low & metadata_baiduII_survival$Sample_ID %in% tcells_patients_high, "C1QC+ CD74+ Macrophages low CD8+ CXCL13+ T cells high","other"))))
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/BaiduII_tumor_macro_tcells_validation_survival_plot.pdf")

km_macrophage_baiduI <- survfit(Surv(Survival_time, Survival_status) ~ Macrophage, data=metadata_baiduII_survival)
summary(km_macrophage_baiduI)
ggsurvplot(km_macrophage_baiduI, data = metadata_baiduII_survival,conf.int = T, pval = T)

km_tcells_baiduI <- survfit(Surv(Survival_time, Survival_status) ~ Tcells, data=metadata_baiduII_survival)
summary(km_tcells_baiduI)
ggsurvplot(km_tcells_baiduI, data = metadata_baiduII_survival,conf.int = T, pval = T)

km_total_baiduI <- survfit(Surv(Survival_time, Survival_status) ~ Total, data=metadata_baiduII_survival)
summary(km_total_baiduI)
ggsurvplot(km_total_baiduI, data = metadata_baiduII_survival,conf.int = T, pval = T)
dev.off()

#### Correlation plot
gsva.res <- as.data.frame(gsva.res)
gsva.res$sample <- rownames(gsva.res)
p <- ggplot(gsva.res,
            aes(x = `C1QC+ CD74+ Macrophages`, y= `CD8+ CXCL13+ T cells`,color =sample )) +
  geom_point() + ylab("CD8+ CXCL13+ T cells") + xlab("C1QC+ CD74+ Macrophages") +theme_classic()
p



##### 5. compare in clinical info #####
## clinical info
gsva.res$patient <- as.numeric(rownames(gsva.res))
metadata_baiduII_merge <- merge(gsva.res,metadata_BaiduII , by.x = "patient", by.y = "Tumor_Sample_Barcode")
## keep
ggplot(metadata_baiduII_merge, aes(x = Grade, y = `C1QC+ CD74+ Macrophages`)) + geom_boxplot()
metadata_baiduII_merge$Grade <- factor(metadata_baiduII_merge$Grade, c("G1","G2","G3"))
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/BaiduIi_clinical_info.pdf")
ggboxplot(metadata_baiduII_merge ,x = "Grade", y = "C1QC+ CD74+ Macrophages" ,fill = "Grade") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("G1","G2"), c("G1","G3"),c("G2","G3")))
dev.off()
## bad
ggplot(metadata_baiduII_merge, aes(x = `TNM stage`, y = `C1QC+ CD74+ Macrophages`)) + geom_boxplot()
## bad
ggplot(metadata_baiduII_merge, aes(x = Metastasis_recurrence_status, y = `C1QC+ CD74+ Macrophages`)) + geom_boxplot()





##### 15. TCGA validation ######
## ESCC
count <- as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/TCGA_counts.csv"))
count <- count[!is.na(count$X),]
count <- count[!duplicated(count$X),]

count %>% distinct(X,.keep_all = T) %>% as.tibble() %>% column_to_rownames(var = "X")
rownames(count) <- count$X
count <- count[,colnames(count)[2:186]]
colnames(count) <- gsub("\\.", "-",colnames(count))
metadata <- as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/TCGA_clinical_info.csv"))
count <- count[,metadata$Sample.ID]
dim(count )

tpm <- read.csv("/Users/lexiiiii/new_onedrive/OneDrive/Work/phD/phd_project/TME_gender/results/TCGA_individual/TCGA_tpm.csv")
## Get CD8+ CXCL13+ Tex cells markers 
cd8_cxcl13_markers <- read.csv("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/res_0.5_cosg_cluster_marker.csv")
cd8_cxcl13_markers <- c(cd8_cxcl13_markers$X11,"CD8A","CXCL13","CD3E")
## Get C1QC+ CD74+ Macrophages markers
c1qc_TAM_markers <- read.csv("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/macrophage/res_2_mcells_cosg_cluster_marker.csv")
c1qc_TAM_markers <- c(c1qc_TAM_markers$X4,"CD74","LYZ")
## GSVA 
gsva_marker_list <- list(c1qc_TAM_markers,cd8_cxcl13_markers )
names(gsva_marker_list) <- c("C1QC+ CD74+ Macrophages","CD8+ CXCL13+ T cells")

count_tcga <- as.matrix(count)
gsva.res <- gsva(count_tcga, gsva_marker_list, method = "gsva")
gsva.res <- as.data.frame(t(gsva.res))
## High group vs low group 
### macrophage ##
macro_patients_high <- rownames(gsva.res[gsva.res[,1] >= mean(gsva.res[,1]),])
#macro_patients_high <- as.numeric(macro_patients_high )
macro_patients_low <- rownames(gsva.res[gsva.res[,1] < mean(gsva.res[,1]),])

### Tcells ##
tcells_patients_high <- rownames(gsva.res[gsva.res[,2] >= mean(gsva.res[,2]),])
#tcells_patients_high <- as.numeric(tcell_patients_high )
tcells_patients_low <- rownames(gsva.res[gsva.res[,2] < mean(gsva.res[,2]),])


### Survival plot 
metadata_TCGA_survival <- metadata
metadata_TCGA_survival$macro_group <- ifelse(metadata_TCGA_survival$Sample.ID %in% macro_patients_high,"C1QC+ CD74+ Macrophages high","C1QC+ CD74+ Macrophages low")
metadata_TCGA_survival$tcells_group <- ifelse(metadata_TCGA_survival$Sample.ID %in% tcells_patients_high, "CD8+ CXCL13+ T cells high","CD8+ CXCL13+ T cells low")
metadata_TCGA_survival$Total <- ifelse(metadata_TCGA_survival$Sample.ID %in% macro_patients_high & metadata_TCGA_survival$Sample.ID %in% tcells_patients_high, "C1QC+ CD74+ Macrophages high CD8+ CXCL13+ T cells high",
                                          ifelse(metadata_TCGA_survival$Sample.ID %in% macro_patients_high & metadata_TCGA_survival$Sample.ID %in% tcells_patients_low, "C1QC+ CD74+ Macrophages high CD8+ CXCL13+ T cells low",
                                                 ifelse(metadata_TCGA_survival$Sample.ID %in% macro_patients_low & metadata_TCGA_survival$Sample.ID %in% tcells_patients_low, "C1QC+ CD74+ Macrophages low CD8+ CXCL13+ T cells low",
                                                        ifelse(metadata_TCGA_survival$Sample.ID %in% macro_patients_low & metadata_TCGA_survival$Sample.ID %in% tcells_patients_high, "C1QC+ CD74+ Macrophages low CD8+ CXCL13+ T cells high","other"))))


pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/TCGA_tumor_macro_tcells_validation_survival_plot.pdf")
metadata_TCGA_survival$survival_status <- ifelse(metadata_TCGA_survival$Overall.Survival.Status == "1:DECEASED",1,0)
macro_km_TCGA_fit <- survfit(Surv(`Overall.Survival..Months.`, survival_status) ~ macro_group, data=metadata_TCGA_survival)
ggsurvplot(macro_km_TCGA_fit, linetype = "strata", pval = T,
           palette = c("#E7B800", "#2E9FDF"),  data = metadata_TCGA_survival[,c("Overall.Survival..Months.","survival_status","macro_group")], risk.table = F)

tcells_km_TCGA_fit <- survfit(Surv(`Overall.Survival..Months.`, survival_status) ~ tcells_group, data=metadata_TCGA_survival)
ggsurvplot(tcells_km_TCGA_fit, linetype = "strata", pval = T,
           palette = c("#E7B800", "#2E9FDF"),  data = metadata_TCGA_survival[,c("Overall.Survival..Months.","survival_status","tcells_group")], risk.table = F)

total_km_TCGA_fit <- survfit(Surv(`Overall.Survival..Months.`, survival_status) ~ Total, data=metadata_TCGA_survival)
ggsurvplot(total_km_TCGA_fit, linetype = "strata", pval = T,
           #palette = c("#E7B800", "#2E9FDF"), 
           data = metadata_TCGA_survival[,c("Overall.Survival..Months.","survival_status","Total")], risk.table = F)


dev.off()

## clinical info
gsva.res$patient <- rownames(gsva.res)
metadata_TCGA_merge <- merge(gsva.res, metadata_TCGA_survival, by.x = "patient", by.y = "Sample.ID")
## bad
metadata_TCGA_merge_keep <- metadata_TCGA_merge[!is.na(metadata_TCGA_merge$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code ),]
metadata_TCGA_merge_keep <- metadata_TCGA_merge[metadata_TCGA_merge$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code %in% c("M0","M1"),]
ggplot(metadata_TCGA_merge_keep, aes(x = American.Joint.Committee.on.Cancer.Metastasis.Stage.Code, y = `C1QC+ CD74+ Macrophages`)) + geom_boxplot()

## Lymph node
metadata_TCGA_merge_keep <- metadata_TCGA_merge[ metadata_TCGA_merge$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code != "NX",]
metadata_TCGA_merge_keep <- metadata_TCGA_merge_keep[ ! is.na(metadata_TCGA_merge$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code),]
metadata_TCGA_merge_keep$lymph_node_stage <- ifelse(metadata_TCGA_merge_keep$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code == "N3","N3","N0,N1,N2")
ggplot(metadata_TCGA_merge_keep, aes(x = Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code, y = `C1QC+ CD74+ Macrophages`)) + geom_boxplot()
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/TCGA_clinical_info.pdf")

ggboxplot(metadata_TCGA_merge_keep ,x = "lymph_node_stage", y = "C1QC+ CD74+ Macrophages" ,
          fill = "lymph_node_stage") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("N0,N1,N2",'N3')))


## bad
ggplot(metadata_TCGA_merge, aes(x = `C1QC+ CD74+ Macrophages`, y = Percent.Lymphocyte.Infiltration)) + geom_point()
## bad
ggplot(metadata_TCGA_merge, aes(x = `C1QC+ CD74+ Macrophages`, y = TNM.Stage)) + geom_boxplot()
## keep 
metadata_TCGA_merge_keep <- metadata_TCGA_merge[!is.na(metadata_TCGA_merge$American.Joint.Committee.on.Cancer.Tumor.Stage.Code ),]
metadata_TCGA_merge_keep$American.Joint.Committee.on.Cancer.Tumor.Stage.Code <- factor(metadata_TCGA_merge_keep$American.Joint.Committee.on.Cancer.Tumor.Stage.Code, c("T1","T2","T3","T4"))
#pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/TCGA_clinical_info.pdf")
ggboxplot(metadata_TCGA_merge_keep ,x = "American.Joint.Committee.on.Cancer.Tumor.Stage.Code", y = "C1QC+ CD74+ Macrophages" ,fill = "American.Joint.Committee.on.Cancer.Tumor.Stage.Code") + scale_fill_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(comparisons = list(c("T1","T2"), c("T2","T3"),c("T3","T4")))
dev.off()



####### 16. cellchat Overlap pathway #####
### BCC
BCC_cellchat <- readRDS("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_only8/cellchat_all/BCC_cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(BCC_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8_ex_T_cells" & c1qc_cd8_scores$target == "C1QC+ TAM",]
BCC_pathway <- unique(c1qc_cd8_scores$pathway_name)

### UVM
UVM_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM/cellchat_all/UVM_cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(UVM_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ CXCL13+ Tex" & c1qc_cd8_scores$target == "C1QC+ TAM",]
UVM_pathway <- unique(c1qc_cd8_scores$pathway_name)

### ESCC_SiT
sit_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(sit_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ Effctor T" & c1qc_cd8_scores$target == "C1QC+ TAM",]
sit_pathway <- unique(c1qc_cd8_scores$pathway_name)

### ESCC GSE203067
GSE_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/validation/ESCC/cellchat_only_LN/only_LN_ESCC_cellchat.rds")
c1qc_cd8_scores <- subsetCommunication(GSE_cellchat)
c1qc_cd8_scores <- c1qc_cd8_scores[c1qc_cd8_scores$source == "CD8+ CXCL13+ Tcells" & c1qc_cd8_scores$target == "C1QC+ TAMs",]
GSE_pathway <- unique(c1qc_cd8_scores$pathway_name)

list_pathway <- list("BCC" =BCC_pathway, "UVM" = UVM_pathway, "ESCC_SiT" = sit_pathway, "ESCC_GSE203067" = GSE_pathway)

venn(list_pathway)
venn.plot <- venn.diagram(
  x = list_pathway,
  filename = "~/onedrive/Work/phD/phd_project/SiT/results/validation/four_groups_venn.PDF",
  col = "black",
  #lty = "dotted",
  fill = c("#300A5BFF", "#C43C4EFF", "#FCA007FF","#ED6925FF"),
  alpha = 0.9,
  #label.col = c("darkred", "white", "darkblue", "white",
   #             "white", "white", "darkgreen"),
  cex = 2.5,#内标签的字体大小
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "outer",#设置标签在圆外面
  #cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2,#外标签的字体大小
  cat.fontfamily = "serif",
  cat.dist = c(0.05, 0.05, 0.05,0.05),#相对圆圈的位置
  cat.pos = c(0,-90,90,180)  #相对12点方向旋转的角度
)
write.table(BCC_pathway, "~/onedrive/Work/phD/phd_project/SiT/results/validation/BCC_pathway.txt",row.names = F, quote = F)
write.table(UVM_pathway, "~/onedrive/Work/phD/phd_project/SiT/results/validation/UVM_pathway.txt",row.names = F, quote = F)
write.table(sit_pathway, "~/onedrive/Work/phD/phd_project/SiT/results/validation/sit_pathway.txt",row.names = F, quote = F)
write.table(GSE_pathway, "~/onedrive/Work/phD/phd_project/SiT/results/validation/GSE_pathway.txt",row.names = F, quote = F)

all_overlap <- intersect(intersect(intersect(BCC_pathway,UVM_pathway),sit_pathway),GSE_pathway)
setdiff(sit_pathway,all_overlap)
setdiff(BCC_pathway,all_overlap)
setdiff(intersect(sit_pathway,BCC_pathway),sit_pathway)