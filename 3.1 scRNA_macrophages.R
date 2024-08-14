######## 0. library import #########
library(Seurat)
library(stringr)
library(harmony)
library(MAST)
library(SingleR)
library(ggplot2)
library(COSG)
library(harmony)
library(DoubletFinder)
library(monocle3)
library(msigdbr)
library(tidyverse)
library(patchwork)
library(FlexDotPlot)
library(GSVA)
library(pheatmap)
library(viridis)
######## 1. Data import ##########
sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")

######## 2. Extract Macrophage only ######
mcells_only_sce_anno_ori <- subset(sce_anno, subset = seurat_clusters == "9")

before_mcells_only_sce_anno <- subset(mcells_only_sce_anno_ori, subset = orig.ident == "before")
after_mcells_only_sce_anno <- subset(mcells_only_sce_anno_ori, subset = orig.ident == "after")
after_ln_mcells_only_sce_anno <- subset(mcells_only_sce_anno_ori, subset = orig.ident == "after_ln")

mcells_only_sce_anno <- SCTransform(mcells_only_sce_anno_ori) %>% RunPCA() 
mcells_only_sce_anno <- RunHarmony(mcells_only_sce_anno, group.by.vars = "orig.ident",
                                   assay.use = "SCT", max.iter.harmony = 10)  
## 降维聚类
pc.num=1:30
mcells_only_sce_anno <- RunUMAP(mcells_only_sce_anno, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=1)   

saveRDS(mcells_only_sce_anno,"~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/mcells_only_sce_anno.rds")
mcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/mcells_only_sce_anno.rds")
######## 3. TAM annotation and extract TAM  (Fig 2C) #######
mcells_only_sce_anno$tam_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(1,2,4),"C1QC+ TAM",
                                         ifelse(mcells_only_sce_anno$seurat_clusters %in% c(0,3,5,6),"FCN1+ TAM","deleted"))

mcells_only_sce_anno$tam_label <- ifelse(mcells_only_sce_anno$seurat_clusters %in% c(1,2,4),"C1QC+ TAM",
                                         ifelse(mcells_only_sce_anno$seurat_clusters %in% c(0,5,6),"FCN1+ TAM",
                                                ifelse(mcells_only_sce_anno$seurat_clusters %in% c(3),"Other TAM","deleted")))


pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/macrophages_all_cluster_plot.pdf",height = 10)
DimPlot(mcells_only_sce_anno, group.by = "seurat_clusters",reduction = "harmony" )
VlnPlot(mcells_only_sce_anno, features =c("C1QC", "FCN1", "CCL5", "SIRPA","LYZ", "CD8A","CD3E","CD74"), group.by = "seurat_clusters")

VlnPlot(mcells_only_sce_anno, features =c("C1QC", "FCN1", "CCL5", "SIRPA","LYZ", "CD8A","CD3E","CD74"), group.by = "tam_label")
DimPlot(mcells_only_sce_anno, group.by = "tam_label",reduction = "harmony" )
dev.off()

saveRDS(mcells_only_sce_anno,"~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/mcells_only_sce_anno_labeled.rds")
mcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/mcells_only_sce_anno_labeled.rds")
#### Extract only TAM

TAM_only_sce_anno <- subset(mcells_only_sce_anno, subset = tam_label %in% c("C1QC+ TAM", "FCN1+ TAM" ))

saveRDS(TAM_only_sce_anno,"~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_only_sce_anno.rds")
TAM_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_only_sce_anno.rds")


########## 4. Marker visualization #########
### dimplot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_dimplot.pdf",width = 8)

DimPlot(TAM_only_sce_anno,  group.by = "tam_label", label = F,reduction = "harmony") 
DimPlot(TAM_only_sce_anno, reduction = "harmony", group.by = "orig.ident", label = F) 
DimPlot(TAM_only_sce_anno, reduction = "harmony", group.by = "seurat_clusters", label = F) 

dev.off()
### marker featureplot
TAM_markers <- c("CD68","LYZ","CSF1R", "MARCO","CD40", "CXCL2", 
                 "CD86","CD80","IL1A","IL6", "IL12", "IL23","FCGR3A", "IFNG", "TNF",
                 "INOS","CD163","MRC1" ,"IL10", "CCL5", "CD81","CIITA","TGFB1"
)
cluster_markers <- c("C1QA", "C1QC","C1QB", "CD74", "FCN1", "APOBEC3A", "ECE1")
m1_like_markers <- c("CXCL9", "CXCL10", "CXCL11",  "TNF","IL12" )
m2_like_markers <- c("ARG1", "YM1", "FIZZ", "CCL2", "CCL17", "CCL22")
total_markers <- c(cluster_markers,m1_like_markers, m2_like_markers)
TAM_only_sce_anno@active.assay <- "RNA"

##### Feature plot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_markers_featureplot.pdf", width = 15, height = 3)
FeaturePlot(TAM_only_sce_anno,reduction = "harmony", features = c("C1QC", "FCN1", "SIRPA","CD74"), pt.size = 0.6,ncol = 4 )+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off( )

### dotplot 
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_markers_dotplot.pdf", width = 8, height = 5)
DotPlot(TAM_only_sce_anno,features = cluster_markers, group.by = "tam_label") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
DotPlot(TAM_only_sce_anno,features = m1_like_markers, group.by = "tam_label") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
DotPlot(TAM_only_sce_anno,features = m2_like_markers, group.by = "tam_label") + 
  scale_color_viridis(begin = 0.2,end = 0.9,option = "E") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5))
dev.off()

## Vlnplot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_markers_vlnplot.pdf")
VlnPlot(TAM_only_sce_anno, features = c("C1QC", "FCN1", "SIRPA","LYZ","CD74"), group.by = "tam_label", cols = c("C1QC+ TAM" = "#440154FF", "FCN1+ TAM" = "#22A884FF")) 
dev.off()
########## 5. GSVA hallmarks ########
####### all
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_genesets <- subset(hallmark_genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
hallmark_genesets <- split(hallmark_genesets$gene_symbol, hallmark_genesets$gs_name)
## total
expr <- AverageExpression(TAM_only_sce_anno, assays = "RNA", slot = "data",group.by = "tam_label")
expr <- as.matrix(expr$RNA)
gsva.res <- gsva(expr, hallmark_genesets, method = "gsva")
rownames(gsva.res) <- gsub("HALLMARK_","",rownames(gsva.res))
colnames(gsva.res) <- c("C1QC+ TAM", "FCN1+ TAM")

annotation_row <-  as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/hallmark.csv",col.names = F))
colnames(annotation_row) <- "Pathway"

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_heatmap_hallmark.pdf",width = 10, height = 10)
pheatmap(gsva.res[rownames(annotation_row),],show_colnames = T, angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()
########## 6. GSVA inflammation (all cells) ########

inflammation_factor <- c("TNF", "IL1A", "IL6", "CXCL8", "IL12RB1")
anti_inflammation_factor <- c("IL10", "TGFB1","IL4","CSF1")
chemokines <- c("CCL2", "CCL5", "CXCL12")
TAM_related_list <- list(inflammation_factor  = inflammation_factor ,
                         anti_inflammation_factor =anti_inflammation_factor,
                         chemokines =chemokines)
gsva.res_other <- gsva( as.matrix(TAM_only_sce_anno@assays$SCT@counts), TAM_related_list, method = "gsva")

### 1) total
gsva.res_other_boxplot <- as.data.frame(t(gsva.res_other))
gsva.res_other_boxplot$id <- rownames(gsva.res_other_boxplot)
gsva.res_other_boxplot <- merge(gsva.res_other_boxplot, TAM_only_sce_anno@meta.data, by.x = "id", by.y = "barcode")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_RNA_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
dev.off()

### 2) total in each sample
gsva.res_other_boxplot <- as.data.frame(t(gsva.res_other))
gsva.res_other_boxplot$id <- rownames(gsva.res_other_boxplot)
gsva.res_other_boxplot <- merge(gsva.res_other_boxplot, TAM_only_sce_anno@meta.data, by.x = "id", by.y = "barcode")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_before_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_after_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_after_ln_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()


### 2) C1QC+ Before vs after vs after_ln
TAM_labled_c1qc <- subset(TAM_only_sce_anno,subset = tam_label == "C1QC+ TAM")
expr <- AverageExpression(TAM_labled_c1qc ,group.by = c("orig.ident"), assays = "SCT", slot = "data")
expr <- as.matrix(expr$SCT)

gsva.res_TAM_c1qc <- gsva( TAM_labled_c1qc@assays$SCT@counts, TAM_related_list, method = "gsva")
gsva.res_TAM_c1qc_boxplot <- as.data.frame(t(gsva.res_TAM_c1qc))
gsva.res_TAM_c1qc_boxplot$barcode <- rownames(gsva.res_TAM_c1qc_boxplot)
c1qc_boxplot_df <- merge(gsva.res_TAM_c1qc_boxplot, TAM_only_sce_anno@meta.data[,c("barcode","orig.ident")], by = "barcode")
c1qc_boxplot_df$orig.ident <- factor(c1qc_boxplot_df$orig.ident, c("before", "after", "after_ln"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_C1QC_samples_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "inflammation_factor", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln"))) 

ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "anti_inflammation_factor", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln")))

ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "chemokines", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln")))

dev.off()


########## 7. GSVA inflammation (groupby sample ) ########
inflammation_factor <- c("TNF", "IL1A", "IL6", "CXCL8", "IL12RB1")
anti_inflammation_factor <- c("IL10", "TGFB1","IL4","CSF1")
chemokines <- c("CCL2", "CCL5", "CXCL12")
TAM_related_list <- list(inflammation_factor  = inflammation_factor ,
                         anti_inflammation_factor =anti_inflammation_factor,
                         chemokines =chemokines)
TAM_only_sce_anno_c1qc <- subset(TAM_only_sce_anno,subset = tam_label == "C1QC+ TAM")
expr <- AverageExpression(TAM_only_sce_anno_c1qc,group.by = c("orig.ident"), assays = "SCT", slot = "data")
expr <- as.matrix(expr$SCT)

gsva.res_TAM_c1qc <- gsva( expr, TAM_related_list, method = "gsva")
gsva.res_other <- gsva( as.matrix(TAM_only_sce_anno@assays$SCT@counts), TAM_related_list, method = "gsva")

### 1) total
gsva.res_other_boxplot <- as.data.frame(t(gsva.res_other))
gsva.res_other_boxplot$id <- rownames(gsva.res_other_boxplot)
gsva.res_other_boxplot <- merge(gsva.res_other_boxplot, TAM_only_sce_anno@meta.data, by.x = "id", by.y = "barcode")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")
ggboxplot(gsva.res_other_boxplot ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test")

dev.off()

### 2) total in each sample
gsva.res_other_boxplot <- as.data.frame(t(gsva.res_other))
gsva.res_other_boxplot$id <- rownames(gsva.res_other_boxplot)
gsva.res_other_boxplot <- merge(gsva.res_other_boxplot, TAM_only_sce_anno@meta.data, by.x = "id", by.y = "barcode")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_before_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "before",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_after_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_after_ln_SCT_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "anti_inflammation_factor", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 

ggboxplot(gsva.res_other_boxplot[gsva.res_other_boxplot$orig.ident == "after_ln",] ,x = "tam_label", y = "chemokines", color = "tam_label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test") 
dev.off()


### 2) C1QC+ Before vs after vs after_ln
TAM_labled_c1qc <- subset(TAM_only_sce_anno,subset = tam_label == "C1QC+ TAM")
expr <- AverageExpression(TAM_labled_c1qc ,group.by = c("orig.ident"), assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

gsva.res_TAM_c1qc <- gsva( TAM_labled_c1qc@assays$SCT@counts, TAM_related_list, method = "gsva")
gsva.res_TAM_c1qc_boxplot <- as.data.frame(t(gsva.res_TAM_c1qc))
gsva.res_TAM_c1qc_boxplot$barcode <- rownames(gsva.res_TAM_c1qc_boxplot)
c1qc_boxplot_df <- merge(gsva.res_TAM_c1qc_boxplot, TAM_only_sce_anno@meta.data[,c("barcode","orig.ident")], by = "barcode")
c1qc_boxplot_df$orig.ident <- factor(c1qc_boxplot_df$orig.ident, c("before", "after", "after_ln"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_C1QC_samples_boxplot_inflammation_scores.pdf",width = 10, height = 10)
ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "inflammation_factor", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln"))) 

ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "anti_inflammation_factor", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln")))

ggboxplot(c1qc_boxplot_df ,x = "orig.ident", y = "chemokines", color = "orig.ident") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("before","after"),c("before","after_ln")))

dev.off()


########## 8. Add label into original sc ########
TAM_label <- mcells_only_sce_anno@meta.data[,"tam_label"]
names(TAM_label ) <-  rownames(mcells_only_sce_anno@meta.data)

sce_anno_tam_labeled <- AddMetaData(
  object = sce_anno,
  metadata = as.data.frame(TAM_label)
)

head(sce_anno_tam_labeled@meta.data)
sce_anno_tam_labeled$celltype_combined <- ifelse(sce_anno_tam_labeled$TAM_label %in% c("C1QC+ TAM","deleted","FCN1+ TAM","Other TAM"),  sce_anno_tam_labeled$TAM_label,sce_anno_tam_labeled$celltype)

table(sce_anno_tam_labeled$celltype_combined)
table(sce_anno_tam_labeled$celltype)


sce_anno_tam_labeled <- subset(sce_anno_tam_labeled, subset = celltype_combined != "deleted")
saveRDS(sce_anno_tam_labeled,"~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sce_anno_tam_labeled.rds")
sce_anno_tam_labeled <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sce_anno_tam_labeled.rds")
########## 9. Correlation plot  #######
### All celltypes
expr <- AverageExpression(sce_anno_tam_labeled,group.by = "celltype_combined", assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

top1000_genes <- names(tail(sort(apply(expr,1,sd)),1000))
#View(expr[top1000_genes,])
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/correlation_TAM_tcells_celltypes_heatmap.pdf",width = 6, height = 8)

pheatmap(cor(expr[top1000_genes,], method = "spearman"))
dev.off()

#### Only TAM to CD8 Tcells
expr_selected <- expr[,c("C1QC+ TAM", "FCN1+ TAM","CD8+ Effctor T" )]
top1000_genes <- names(tail(sort(apply(expr_selected ,1,sd)),1000))

cor_selected <- as.data.frame(cor(expr_selected[top1000_genes,], method = "spearman"))
cor_selected_df <- data.frame(PCC =cor_selected$`CD8+ Effctor T`, celltype = rownames(cor_selected) )
cor_selected_df <- cor_selected_df[cor_selected_df$celltype != "CD8+ Effctor T",]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/correlation_TAM_tcells_cd8_only_lolipop.pdf",width = 6, height = 4)

ggdotchart(cor_selected_df, y = "PCC", x ="celltype", color = "celltype",dot.size = 6,
           add = "segments") 
dev.off()


########## 10. C1QC tam ratio vs CD8 effector ratio ####
ratio_metadata <- sce_anno@meta.data[,c("celltype", "TAM_label", "orig.ident","celltype_combined")]

cell_ratio_sample_total <- prop.table(table(sce_anno_tam_labeled$celltype_combined,sce_anno_tam_labeled$orig.ident), margin = 2)
#cell_ratio_sample_total <- prop.table(table(ratio_metadata$combined_label,ratio_metadata$orig.ident), margin = 2)

cell_ratio_sample_total

dotplot_df_tam <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                             "CD8_CXCl13_Tcell_ratio" =  cell_ratio_sample_total["CD8+ CXCL13+ Tex",],
                             sample = colnames(cell_ratio_sample_total))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_cd8_ratio_plot.pdf",width = 5, height =4)
ggplot(dotplot_df_tam, aes(x=C1QC_TAM_ratio, y= `CD8_CXCl13_Tcell_ratio`, size = 50,color=sample)) +
  geom_point(alpha=1) + scale_size(range = c(2, 15)) + scale_color_viridis(discrete = T)
dev.off()




########## 11. FCN1 vs C1QC Volcano plot (Fig 2D) ####


########## 12. Sep Cellchat ( C1qc TAM and CD8 tex) (Fig 2J,K) ########

sce_anno_tam_labeled_focus <- subset(sce_anno_tam_labeled, subset = celltype_combined %in% c("C1QC+ TAM", 'CD8+ Effctor T' ))

before_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "before")
after_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "after")
after_ln_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "after_ln")

setwd("cellchat_focus")

### 11.1 create cellchat 
before_cellchat <- createCellChat(before_t_TAM_sce, group.by = "celltype_combined")
after_cellchat <- createCellChat(after_t_TAM_sce, group.by = "celltype_combined")
after_ln_cellchat <- createCellChat(after_ln_t_TAM_sce, group.by = "celltype_combined")

### 11.2 构建细胞通讯网络
## before
before_cellchat@DB <- CellChatDB.human
before_cellchat <- subsetData(before_cellchat) # This step is necessary even if using the whole database
before_cellchat <- identifyOverExpressedGenes(before_cellchat)
before_cellchat <- identifyOverExpressedInteractions(before_cellchat)

before_cellchat <- computeCommunProb(before_cellchat)
before_cellchat <- filterCommunication(before_cellchat, min.cells = 10)
before_cellchat <- computeCommunProbPathway(before_cellchat)
before_cellchat <- aggregateNet(before_cellchat)

saveRDS(before_cellchat, "before_cellchat.rds")
#before_cellchat <= before_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/before_cellchat.rds")
## after
after_cellchat@DB <- CellChatDB.human
after_cellchat <- subsetData(after_cellchat) # This step is necessary even if using the whole database
after_cellchat <- identifyOverExpressedGenes(after_cellchat)
after_cellchat <- identifyOverExpressedInteractions(after_cellchat)

after_cellchat <- computeCommunProb(after_cellchat)
after_cellchat <- filterCommunication(after_cellchat, min.cells = 10)
after_cellchat <- computeCommunProbPathway(after_cellchat)
after_cellchat <- aggregateNet(after_cellchat)

saveRDS(after_cellchat, "after_cellchat.rds")
#after_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/after_cellchat.rds")

## after_ln
after_ln_cellchat@DB <- CellChatDB.human
after_ln_cellchat <- subsetData(after_ln_cellchat) # This step is necessary even if using the whole database
after_ln_cellchat <- identifyOverExpressedGenes(after_ln_cellchat)
after_ln_cellchat <- identifyOverExpressedInteractions(after_ln_cellchat)

after_ln_cellchat <- computeCommunProb(after_ln_cellchat)
after_ln_cellchat <- filterCommunication(after_ln_cellchat, min.cells = 10)
after_ln_cellchat <- computeCommunProbPathway(after_ln_cellchat)
after_ln_cellchat <- aggregateNet(after_ln_cellchat)

saveRDS(after_ln_cellchat, "after_ln_cellchat.rds")
#after_ln_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/after_ln_cellchat.rds")

### 11.3 Merge CellChat object
before_after_cellchat_list <- list(before = before_cellchat, after = after_cellchat)
before_after_merged_cellchat <- mergeCellChat(before_after_cellchat_list, add.names = names(before_after_cellchat_list),
                                              cell.prefix = T)
before_after_ln_cellchat_list <- list(before = before_cellchat, after_ln = after_ln_cellchat)
before_after_ln_merged_cellchat <- mergeCellChat(before_after_ln_cellchat_list, add.names = names(before_after_ln_cellchat_list),
                                                 cell.prefix = T)


saveRDS(before_after_merged_cellchat,"before_after_merged_cellchat.rds")
saveRDS(before_after_ln_merged_cellchat,"before_after_ln_merged_cellchat.rds")

### 11.4 可视化

# （1） 通讯数量与强度对比
gg1 <- compareInteractions(before_after_merged_cellchat, show.legend = F , group = c(1,2))
gg2 <- compareInteractions(before_after_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("before_after_compareInteractions.pdf")
gg1 + gg2
dev.off()

gg1 <- compareInteractions(before_after_ln_merged_cellchat, show.legend = F,group = c(1,2))
gg2 <- compareInteractions(before_after_ln_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("before_after_ln_compareInteractions.pdf")
gg1 + gg2
dev.off()

# (2) 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_groups_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(before_after_merged_cellchat, weight.scale = T, measure = "weight")
netVisual_diffInteraction(before_after_ln_merged_cellchat, weight.scale = T, measure = "weight")
#netVisual_diffInteraction(total_merged_cellchat, weight.scale = T, measure = "weight")

dev.off()


# (3) Heatmap
pdf("before_after_netVisual_heatmap.pdf")
netVisual_heatmap(before_after_merged_cellchat, measure = "weight")
dev.off()

pdf("before_after_ln_netVisual_heatmap.pdf")
netVisual_heatmap(before_after_ln_merged_cellchat, measure = "weight")
dev.off()

#pdf("all_netVisual_heatmap.pdf")
#netVisual_heatmap(total_merged_cellchat, measure = "weight")
#dev.off()

# (4) 特定信号通路的对比
# weight.max <- getMaxWeight(total_cellchat_list , attribute = c("idents","count"))
# par(mfrow = c(1,2), xpd=TRUE)
# pdf("all_netVisual_circle.pdf")
# for (i in 1:length(total_cellchat_list )) {
#   netVisual_circle(total_cellchat_list [[i]]@net$count, weight.scale = T,
#                    label.edge= F, edge.weight.max = weight.max[2], 
#                    edge.width.max = 12, title.name = paste0("Number of interactions - ", 
#                                                             names(total_cellchat_list )[i]))
# }
# dev.off()

pdf("MIF_netVisual_circle.pdf")
for (i in 1:length(total_cellchat_list )) {
  netVisual_aggregate(total_cellchat_list[[i]], layout = "circle", signaling = c("MIF"),
                      label.edge= F, edge.weight.max = weight.max[2], 
                      edge.width.max = 12, signaling.name = paste0("MIF  - ", 
                                                                   names(total_cellchat_list )[i]))
}
dev.off()

# (5) 通路信号强度对比分析

pdf("before_after_rankNet.pdf")
rankNet(before_after_merged_cellchat, mode= "comparison",stacked = T)
dev.off()

pdf("before_after_ln_rankNet.pdf")
rankNet(before_after_ln_merged_cellchat, mode= "comparison",stacked = T)
dev.off()

# (6) 配体-受体对比分析
pdf("after_TAM_cd8_total_netVisual_bubble.pdf",height = 10)
netVisual_bubble(before_after_merged_cellchat,
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in before", max.dataset = 1)
netVisual_bubble(before_after_merged_cellchat, 
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Lower in before", max.dataset = 2)
dev.off()



pdf("after_ln_TAM_cd8_total_netVisual_bubble.pdf",height = 10)

netVisual_bubble(before_after_ln_merged_cellchat,
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in before", max.dataset = 1)
netVisual_bubble(before_after_ln_merged_cellchat, 
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Lower in before", max.dataset = 2)
dev.off()







########## 13. Sep Cellchat ( all TAM and all Tcells)########
tcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcells_only_sce_anno_ori.rds")
sce_anno_tam_labeled_focus <- subset(sce_anno_tam_labeled, subset = celltype_combined %in% c("C1QC+ TAM", "FCN1+ TAM",unique(tcells_only_sce_anno$celltype) ))
sce_anno_tam_labeled_focus <- SCTransform(sce_anno_tam_labeled_focus) %>% RunPCA() 

before_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "before")
after_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "after")
after_ln_t_TAM_sce <- subset(sce_anno_tam_labeled_focus, subset = orig.ident == "after_ln")

setwd("cellchat_tcells_tam")

### 11.1 create cellchat 
before_cellchat <- createCellChat(before_t_TAM_sce, group.by = "celltype_combined")
after_cellchat <- createCellChat(after_t_TAM_sce, group.by = "celltype_combined")
after_ln_cellchat <- createCellChat(after_ln_t_TAM_sce, group.by = "celltype_combined")

### 11.2 构建细胞通讯网络
## before
before_cellchat@DB <- CellChatDB.human
before_cellchat <- subsetData(before_cellchat) # This step is necessary even if using the whole database
before_cellchat <- identifyOverExpressedGenes(before_cellchat)
before_cellchat <- identifyOverExpressedInteractions(before_cellchat)

before_cellchat <- computeCommunProb(before_cellchat)
before_cellchat <- filterCommunication(before_cellchat, min.cells = 10)
before_cellchat <- computeCommunProbPathway(before_cellchat)
before_cellchat <- aggregateNet(before_cellchat)

saveRDS(before_cellchat, "before_cellchat.rds")
#before_cellchat <= before_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/before_cellchat.rds")
## after
after_cellchat@DB <- CellChatDB.human
after_cellchat <- subsetData(after_cellchat) # This step is necessary even if using the whole database
after_cellchat <- identifyOverExpressedGenes(after_cellchat)
after_cellchat <- identifyOverExpressedInteractions(after_cellchat)

after_cellchat <- computeCommunProb(after_cellchat)
after_cellchat <- filterCommunication(after_cellchat, min.cells = 10)
after_cellchat <- computeCommunProbPathway(after_cellchat)
after_cellchat <- aggregateNet(after_cellchat)

saveRDS(after_cellchat, "after_cellchat.rds")
#after_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/after_cellchat.rds")

## after_ln
after_ln_cellchat@DB <- CellChatDB.human
after_ln_cellchat <- subsetData(after_ln_cellchat) # This step is necessary even if using the whole database
after_ln_cellchat <- identifyOverExpressedGenes(after_ln_cellchat)
after_ln_cellchat <- identifyOverExpressedInteractions(after_ln_cellchat)

after_ln_cellchat <- computeCommunProb(after_ln_cellchat)
after_ln_cellchat <- filterCommunication(after_ln_cellchat, min.cells = 10)
after_ln_cellchat <- computeCommunProbPathway(after_ln_cellchat)
after_ln_cellchat <- aggregateNet(after_ln_cellchat)

saveRDS(after_ln_cellchat, "after_ln_cellchat.rds")
#after_ln_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/sep/after_ln_cellchat.rds")

### 11.3 Merge CellChat object
before_after_cellchat_list <- list(before = before_cellchat, after = after_cellchat)
before_after_merged_cellchat <- mergeCellChat(before_after_cellchat_list, add.names = names(before_after_cellchat_list),
                                              cell.prefix = T)
before_after_ln_cellchat_list <- list(before = before_cellchat, after_ln = after_ln_cellchat)
before_after_ln_merged_cellchat <- mergeCellChat(before_after_ln_cellchat_list, add.names = names(before_after_ln_cellchat_list),
                                                 cell.prefix = T)


saveRDS(before_after_merged_cellchat,"before_after_merged_cellchat.rds")
saveRDS(before_after_ln_merged_cellchat,"before_after_ln_merged_cellchat.rds")

### 11.4 可视化

# （1） 通讯数量与强度对比
gg1 <- compareInteractions(before_after_merged_cellchat, show.legend = F , group = c(1,2))
gg2 <- compareInteractions(before_after_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("before_after_compareInteractions.pdf")
gg1 + gg2
dev.off()

gg1 <- compareInteractions(before_after_ln_merged_cellchat, show.legend = F,group = c(1,2))
gg2 <- compareInteractions(before_after_ln_merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("before_after_ln_compareInteractions.pdf")
gg1 + gg2
dev.off()

# (2) 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf("all_groups_netVisual_diffInteraction.pdf")
netVisual_diffInteraction(before_after_merged_cellchat, weight.scale = T, measure = "weight")
netVisual_diffInteraction(before_after_ln_merged_cellchat, weight.scale = T, measure = "weight")
#netVisual_diffInteraction(total_merged_cellchat, weight.scale = T, measure = "weight")

dev.off()


# (3) Heatmap
pdf("before_after_netVisual_heatmap.pdf")
netVisual_heatmap(before_after_merged_cellchat, measure = "weight")
dev.off()

pdf("before_after_ln_netVisual_heatmap.pdf")
netVisual_heatmap(before_after_ln_merged_cellchat, measure = "weight")
dev.off()

#pdf("all_netVisual_heatmap.pdf")
#netVisual_heatmap(total_merged_cellchat, measure = "weight")
#dev.off()

# (4) 特定信号通路的对比
# weight.max <- getMaxWeight(total_cellchat_list , attribute = c("idents","count"))
# par(mfrow = c(1,2), xpd=TRUE)
# pdf("all_netVisual_circle.pdf")
# for (i in 1:length(total_cellchat_list )) {
#   netVisual_circle(total_cellchat_list [[i]]@net$count, weight.scale = T,
#                    label.edge= F, edge.weight.max = weight.max[2], 
#                    edge.width.max = 12, title.name = paste0("Number of interactions - ", 
#                                                             names(total_cellchat_list )[i]))
# }
# dev.off()

pdf("MIF_netVisual_circle.pdf")
for (i in 1:length(total_cellchat_list )) {
  netVisual_aggregate(total_cellchat_list[[i]], layout = "circle", signaling = c("MIF"),
                      label.edge= F, edge.weight.max = weight.max[2], 
                      edge.width.max = 12, signaling.name = paste0("MIF  - ", 
                                                                   names(total_cellchat_list )[i]))
}
dev.off()

# (5) 通路信号强度对比分析

pdf("before_after_rankNet.pdf")
rankNet(before_after_merged_cellchat, mode= "comparison",stacked = T)
dev.off()

pdf("before_after_ln_rankNet.pdf")
rankNet(before_after_ln_merged_cellchat, mode= "comparison",stacked = T)
dev.off()

# (6) 配体-受体对比分析
pdf("after_TAM_cd8_total_netVisual_bubble.pdf",height = 10)
netVisual_bubble(before_after_merged_cellchat,sources.use = c(4,5),targets.use = c(1),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in before", max.dataset = 1)
netVisual_bubble(before_after_merged_cellchat, sources.use = c(4,5),targets.use = c(1),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in after", max.dataset = 2)
dev.off()



pdf("after_ln_TAM_cd8_total_netVisual_bubble.pdf",height = 10)

netVisual_bubble(before_after_ln_merged_cellchat,sources.use = c(4,5),targets.use = c(1),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "Higher in before", max.dataset = 1)
netVisual_bubble(before_after_ln_merged_cellchat, sources.use = c(4,5),targets.use = c(1),
                 comparison = c(1,2), angle.x = 45,  remove.isolate = F,
                 title.name = "higher in after LN", max.dataset = 2)
dev.off()







########## 14. MIF-CD74和 TAM  CD8 CXCL13 ratio correlation ########
## use old
before_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/old_res_0.5/sep/before_cellchat.rds")
after_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/old_res_0.5/sep/after_cellchat.rds")
after_ln_cellchat <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/old_res_0.5/sep/after_ln_cellchat.rds")
### before
before_mif_cd74_scores <- subsetCommunication(before_cellchat)
before_mif_cd74_scores <- before_mif_cd74_scores[before_mif_cd74_scores$source == "CD8+ Effctor T" & before_mif_cd74_scores$target == "C1QC+ TAM",]
before_mif_cd74_scores <- before_mif_cd74_scores[before_mif_cd74_scores$interaction_name == "MIF_CD74_CD44","prob"]
# 0.2196616
### after
after_mif_cd74_scores <- subsetCommunication(after_cellchat)
after_mif_cd74_scores <- after_mif_cd74_scores[after_mif_cd74_scores$source == "CD8+ Effctor T" & after_mif_cd74_scores$target == "C1QC+ TAM",]
after_mif_cd74_scores <- after_mif_cd74_scores[after_mif_cd74_scores$interaction_name == "MIF_CD74_CD44","prob"]
# 0.2153871
### after_ln
after_ln_mif_cd74_scores <- subsetCommunication(after_ln_cellchat)
after_ln_mif_cd74_scores <- after_ln_mif_cd74_scores[after_ln_mif_cd74_scores$source == "CD8+ Effctor T" & after_ln_mif_cd74_scores$target == "C1QC+ TAM",]
after_ln_mif_cd74_scores <- after_ln_mif_cd74_scores[after_ln_mif_cd74_scores$interaction_name == "MIF_CD74_CD44","prob"]
# 0.1891068

### add label into total sc

cell_ratio_sample_total <- prop.table(table(sce_anno_tam_labeled$celltype_combined,sce_anno_tam_labeled$orig.ident), margin = 2)

## 1）Total ratio
dotplot_df <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8+ Effctor T",],
                         cd74_prob = c(after_mif_cd74_scores,after_ln_mif_cd74_scores,
                                       before_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
p1 = ggplot(dotplot_df, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(10, 15)) + scale_color_viridis(discrete = T)


## 2）C1QC / FCN1 TAM ratio
cell_ratio_sample_TAM <- prop.table(table( TAM_only_sce_anno$tam_label,TAM_only_sce_anno$orig.ident), margin = 2)

dotplot_df_tam <- data.frame("C1QC_TAM_FCN1_TAM_ratio" = cell_ratio_sample_TAM["C1QC+ TAM",],
                             "CD8_Effctor_T" = cell_ratio_sample_total["CD8+ Effctor T",],
                             cd74_prob = c(after_mif_cd74_scores,after_ln_mif_cd74_scores,
                                           before_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
p2 = ggplot(dotplot_df_tam, aes(x=C1QC_TAM_FCN1_TAM_ratio, y=CD8_Effctor_T, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(10, 15)) + scale_color_viridis(discrete = T)

## 3）C1QC / Cd8  TAM ratio
dotplot_df_tam <- data.frame("C1QC_TAM_CD8_Effctor_T_ratio" = cell_ratio_sample_total["C1QC+ TAM",] / cell_ratio_sample_total["CD8+ Effctor T",],
                             cd74_prob = c(after_mif_cd74_scores,after_ln_mif_cd74_scores,
                                           before_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))
p3 <- ggplot(dotplot_df_tam, aes(x=C1QC_TAM_CD8_Effctor_T_ratio, y=cd74_prob, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(10, 15)) + scale_color_viridis(discrete = T)

## 4) CD74 exp vs CD8 effector exp
### only CD8 effector
cd8_tex_only <- subset(sce_anno_tam_labeled_focus, subset = celltype_combined == "CD8+ Effctor T")
expr_cd8 <- AverageExpression(cd8_tex_only, assays = "SCT", slot = "data",group.by = "orig.ident")
expr_cd8 <- expr_cd8$SCT["MIF",]

### only C1QC TAM
c1qc_only <- subset(sce_anno_tam_labeled_focus, subset = celltype_combined == "C1QC+ TAM")
expr_c1qc_only <- AverageExpression(c1qc_only, assays = "SCT", slot = "data",group.by = "orig.ident")
expr_c1qc_only <- expr_c1qc_only$SCT["CD74",]

### plot
dotplot_mif_cd74 <- data.frame("MIF_expression" = expr_cd8, "CD74_expression" = expr_c1qc_only,
                               cd74_prob = c(after_mif_cd74_scores,after_ln_mif_cd74_scores,
                                             before_mif_cd74_scores),sample = names(expr_cd8))

p4 <- ggplot(dotplot_mif_cd74, aes(x=MIF_expression, y=CD74_expression, size=cd74_prob, color=sample)) +
  geom_point(alpha=0.5) + scale_size(range = c(10, 15)) + scale_color_viridis(discrete = T)




pdf("scRNA_ratio_MIF_CD74_CD44_pathway_point_plot.pdf",width = 15)
p1+p2+p3+p4
dev.off()
###  5) CD74 exp vs CD8 effector ratio vs C1QC TAM ratio
dotplot_5 <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T" = cell_ratio_sample_total["CD8+ Effctor T",],
                         cd74_expression = expr_c1qc_only,sample = colnames(cell_ratio_sample_total))
p5 <- ggplot(dotplot_5, aes(x=C1QC_TAM_ratio, y=CD8_Effctor_T, size=cd74_expression, color=sample)) +
  geom_point() + scale_size(range = c(8, 15)) + scale_color_viridis(discrete = T) + theme_classic() + theme(legend.position = "top")
pdf("scRNA_ratio_cd74_expression_point_plot.pdf",width= 8) 
p5
dev.off()


### 3D
library(plotly)
dotplot_3d <- data.frame("C1QC_TAM_ratio" = cell_ratio_sample_total["C1QC+ TAM",],
                         "CD8_Effctor_T_ratio" = cell_ratio_sample_total["CD8+ Effctor T",],
                         "CD74_Expression" =  expr_c1qc_only,
                          cd74_prob = c(after_mif_cd74_scores,after_ln_mif_cd74_scores,
                                       before_mif_cd74_scores),sample = colnames(cell_ratio_sample_total))

pdf("scRNA_total_ratio_3d_plot.pdf",width = 15)
p= plot_ly(dotplot_3d, x= ~C1QC_TAM_ratio, y= ~CD8_Effctor_T_ratio, 
            z = ~CD74_Expression, color = ~cd74_prob, marker = list(
             size = 20), name = ~sample) %>% add_markers()
p
dev.off()
