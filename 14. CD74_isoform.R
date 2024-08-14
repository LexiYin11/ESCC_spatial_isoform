########## 0. Library import and functions #####
library(SPATA2)
library(Seurat)
library(ggsci)
library(viridis)
library(msigdbr)
library(GSVA)
library(matrixStats)
library(ggstatsplot)
library(ggpubr)
scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)

}
########## 1. Data import #######
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_added_isoform.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_added_isoform.rds")
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_added_isoform.rds")

after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_ratio_labeled.rds")
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_ratio_labeled.rds")

########## 2. Check CD74 isoform distribution and expression in three samples (Fig 4A,H, Fig S11A) ######
isoforms <- c("CD74..ENST00000353334","CD74..ENST00000009530")
DefaultAssay(after_ln) <- "MULTI"
DefaultAssay(after) <- "MULTI"
DefaultAssay(before) <- "MULTI"
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cd74_isoforms.pdf")
SpatialFeaturePlot(before, features = isoforms , pt.size.factor=3, alpha = c(1, 1))
SpatialFeaturePlot(after, features = isoforms , pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 1) 
SpatialFeaturePlot(after_ln, features = isoforms , pt.size.factor=1.3, alpha = c(1, 1))
dev.off()

########## 3. Extract CD74 isoforms and add into object #####
#### before
isoform1 <- before@assays$MULTI@counts[isoforms[1],]
isoform2 <- before@assays$MULTI@counts[isoforms[2],]
dim(isoform1)
ratio1 <- isoform1 / isoform2
ratio2 <- isoform2 / isoform1
ratio <- rbind(ratio1,ratio2)
dim(ratio)
before[["RATIO"]] <- CreateAssayObject(counts = ratio)
#before <- NormalizeData(object = before, assay = "RATIO")
#before <- ScaleData(before, assay = "RATIO")

DefaultAssay(object = before) <- "RATIO"
pdf("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/cd74_ratio.pdf")
SpatialFeaturePlot(before, features = "ratio1", pt.size.factor=1.8, alpha = c(1, 1),max.cutoff = 15 )
SpatialFeaturePlot(before, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 1)
SpatialFeaturePlot(before, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 4)
SpatialFeaturePlot(before, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 7)

dev.off()

#### after
isoform1 <- after@assays$MULTI@counts[isoforms[1],]
isoform2 <- after@assays$MULTI@counts[isoforms[2],]

ratio1 <- isoform1 / isoform2
ratio2 <- isoform2 / isoform1
summary(ratio1)
ratio <- rbind(ratio1,ratio2)
dim(ratio)
after[["RATIO"]] <- CreateAssayObject(counts = ratio)

DefaultAssay(object = after) <- "RATIO"
pdf("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/cutoff_0.75_ratio_cd74_ratio.pdf")
SpatialFeaturePlot(after, features = "ratio1", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 9)
SpatialFeaturePlot(after, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 1)
SpatialFeaturePlot(after, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 4)
SpatialFeaturePlot(after, features = "ratio2", pt.size.factor=1.8, alpha = c(1, 1),min.cutoff = 7)

dev.off()
#### after_ln
isoform1 <- after_ln@assays$MULTI@counts[isoforms[1],]
isoform2 <- after_ln@assays$MULTI@counts[isoforms[2],]

dim(isoform1)
ratio1 <- isoform1 / isoform2
ratio2 <- isoform2 / isoform1
ratio <- rbind(ratio1,ratio2)
dim(ratio)

after_ln[["RATIO"]] <- CreateAssayObject(counts = ratio)


DefaultAssay(object = after_ln) <- "RATIO"
pdf("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/cutoff_2_cd74_ratio.pdf")
SpatialFeaturePlot(after_ln, features = "ratio1", pt.size.factor=1.3, alpha = c(1, 1),max.cutoff = 15 )
SpatialFeaturePlot(after_ln, features = "ratio2", pt.size.factor=1.3, alpha = c(1, 1),min.cutoff = 1)
SpatialFeaturePlot(after_ln, features = "ratio2", pt.size.factor=1.3, alpha = c(1, 1),min.cutoff = 4)
SpatialFeaturePlot(after_ln, features = "ratio2", pt.size.factor=1.3, alpha = c(1, 1),min.cutoff = 7)

dev.off()



########## 4. Group CD74 ratio of each cell  #####
##### 4.1 before
counts <- t(as.data.frame(before@assays$RATIO@counts))
counts <- na.omit(counts[,"ratio1"])
length(counts)
counts_1 <- counts[counts <= quantile(counts, 0.25) | counts == "Inf"]
length(counts_1)
counts_2 <- counts[counts > quantile(counts, 0.25) & counts <= 5]
length(counts_2)
counts_3 <- counts[counts > 5 & counts != "Inf"]
length(counts_3)

cells_ratio_1 <- names(counts_1)
cells_ratio_2 <- names(counts_2)
cells_ratio_3 <- names(counts_3)
##### 4.2 after
after_counts <- t(as.data.frame(after@assays$RATIO@counts))
after_counts <- na.omit(after_counts[,"ratio1"])
length(after_counts)
after_counts <- ifelse(after_counts == Inf, 0, after_counts)
after_counts_no_0 <- after_counts[after_counts != 0]
quantile(after_counts_no_0)

after_counts_1 <- after_counts[after_counts <= quantile(after_counts, 0.25) | after_counts == "Inf"]
length(after_counts_1)
after_counts_2 <- after_counts[after_counts > quantile(after_counts, 0.25)  & after_counts <= quantile(after_counts, 0.75)]
length(after_counts_2)
after_counts_3 <- after_counts[after_counts > quantile(after_counts, 0.25) & after_counts != "Inf"]
length(after_counts_3)

after_cells_ratio_1 <- names(after_counts_1)
after_cells_ratio_2 <- names(after_counts_2)
after_cells_ratio_3 <- names(after_counts_3)

##### 4.3 after_ln
after_ln_counts <- t(as.data.frame(after_ln@assays$RATIO@counts))
after_ln_counts <- na.omit(after_ln_counts[,"ratio1"])
length(after_ln_counts)

after_ln_counts <- ifelse(after_ln_counts == Inf, 0, counts)
after_ln_counts_no_0 <- counts[after_ln_counts != 0]
quantile(after_ln_counts_no_0)

after_ln_counts_1 <- after_ln_counts[after_ln_counts <= quantile(after_ln_counts, 0.25) | after_ln_counts == "Inf"]
length(after_ln_counts_1)
after_ln_counts_2 <- after_ln_counts[after_ln_counts > quantile(after_ln_counts, 0.25) & after_ln_counts <= quantile(after_counts, 0.75)]
length(after_ln_counts_2)
after_ln_counts_3 <- after_ln_counts[after_ln_counts > quantile(after_counts, 0.75) & after_ln_counts != "Inf"]
length(after_ln_counts_3)

after_ln_cells_ratio_1 <- names(after_ln_counts_1)
after_ln_cells_ratio_2 <- names(after_ln_counts_2)
after_ln_cells_ratio_3 <- names(after_ln_counts_3)

########## 5. Add ratio label (Fig 4B,I,Fig S11B) ####
#### 5.1 before
DefaultAssay(object = before) <- "SCT"
before@meta.data$label <- ifelse(rownames(before@meta.data) %in% cells_ratio_1,"CD74_isoform_ratio_low",
                                ifelse(rownames(before@meta.data) %in% cells_ratio_2,"CD74_isoform_ratio_medium",
                                       ifelse(rownames(before@meta.data) %in% cells_ratio_3,"CD74_isoform_ratio_high","CD74_isoform_ratio_low")))
label_df <- data.frame(Barcode = rownames(before@meta.data), label = before@meta.data$label)

saveRDS(before,"~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_ratio_labeled.rds")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/before_cd74_ratio_dimplot.pdf")
cols <- viridis(3,begin = 0.2,end = 0.8,option = "B")
names(cols) <- c("CD74_isoform_ratio_low", "CD74_isoform_ratio_medium","CD74_isoform_ratio_high")
SpatialDimPlot(before,group.by = "label", pt.size.factor=2.8, alpha = c(1, 1), cols =cols)
#+ scale_fill_viridis(begin = 0,end = 0.4,discrete = F)
dev.off()

#### 5.2 after
DefaultAssay(object = after) <- "SCT"
after@meta.data$label <- ifelse(rownames(after@meta.data) %in% after_cells_ratio_1,"CD74_isoform_ratio_low",
                                   ifelse(rownames(after@meta.data) %in% after_cells_ratio_2,"CD74_isoform_ratio_medium",
                                          ifelse(rownames(after@meta.data) %in% after_cells_ratio_3,"CD74_isoform_ratio_high","CD74_isoform_ratio_low")))
label_df <- data.frame(Barcode = rownames(after@meta.data), label = after@meta.data$label)

saveRDS(after,"~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_ratio_labeled.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_ratio_labeled.rds")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_2_5_after_cd74_ratio_dimplot.pdf")
cols <- c("#1B0C42FF","#CF4446FF", "#FB9A06FF")
names(cols) <- c("CD74_isoform_ratio_low", "CD74_isoform_ratio_medium","CD74_isoform_ratio_high")
SpatialDimPlot(after,group.by = "label", pt.size.factor=1.9, alpha = c(1, 1), cols =cols)
dev.off()


#### 5.3 after_ln
DefaultAssay(object = after_ln) <- "SCT"

after_ln@meta.data$label <- ifelse(rownames(after_ln@meta.data) %in% after_ln_cells_ratio_1,"CD74_isoform_ratio_low",
                                   ifelse(rownames(after_ln@meta.data) %in% after_ln_cells_ratio_2,"CD74_isoform_ratio_medium",
                                          ifelse(rownames(after_ln@meta.data) %in% after_ln_cells_ratio_3,"CD74_isoform_ratio_high","CD74_isoform_ratio_low")))

label_df <- data.frame(Barcode = rownames(after_ln@meta.data), label = after_ln@meta.data$label)
## remove TACCGGTCGTTTCCAT-1
which(label_df$Barcode == "TACCGGTCGTTTCCAT-1")

write.csv(label_df,"~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/ratio_label.csv")

saveRDS(after_ln,"~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_2_5_cd74_ratio_dimplot.pdf")
cols <- c("#1B0C42FF","#CF4446FF", "#FB9A06FF")
names(cols) <- c("low", "medium","high")
SpatialDimPlot(after_ln,group.by = "label", pt.size.factor=1.3, alpha = c(1, 1), cols =cols)
dev.off()

DefaultAssay(after_ln) <- "RATIO"
SpatialFeaturePlot(after_ln, features = "ratio1", pt.size.factor=1.3, alpha = c(1, 1))

########## 6. Compare CD8 CXCL13 expression #####

expr <- AverageExpression(after_ln, assays = "SCT", slot = "data")
expr <- as.matrix(expr$SCT)
expr <- expr[c("CD8A","CXCL13"),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/ratio_label_feature_dotplot.pdf",width = 12, height = 5)

DotPlot(after_ln,features = c("CD8A", "CXCL13"),
        group.by = "label",dot.scale = 7,scale.min = 10,scale = F) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) + scale_color_viridis(option = "A")
dev.off()

########## 7. Create CD8 CXCL13 T cells location feature plot  (Fig 4C,J) ########
## overlap cd8 cxcxl13 cd3e
## 1) after_ln
counts <- t(as.data.frame(after_ln@assays$Spatial@counts))
counts <- counts[,c("CD3E","CD8A","CXCL13")]
dim(counts)

# All > 1
counts_filtered <- counts[rowMins(counts) > 0 ,]
dim(counts_filtered)
after_ln_cd8_cxcl13_spots <- rownames(counts_filtered )
# add label
after_ln@meta.data$cd8_cxcl13_label <- ifelse(rownames(after_ln@meta.data) %in% after_ln_cd8_cxcl13_spots,"CD8_CXCL13_Tcells","Others")
table(after_ln@meta.data$cd8_cxcl13_label)

## 2) after
counts <- t(as.data.frame(after@assays$Spatial@counts))
counts <- counts[,c("CD3E","CD8A","CXCL13")]
dim(counts)
# All > 1
counts_filtered <- counts[rowMins(counts) > 0,]
dim(counts_filtered)
after_cd8_cxcl13_spots <- rownames(counts_filtered )
# add label
after@meta.data$cd8_cxcl13_label <- ifelse(rownames(after@meta.data) %in% after_cd8_cxcl13_spots,"CD8_CXCL13_Tcells","Others")
table(after@meta.data$cd8_cxcl13_label)

## 3) before
counts <- t(as.data.frame(before@assays$Spatial@counts))
counts <- counts[,c("CD8B", "CD3E")]
dim(counts)
# All > 1
counts_filtered <- counts[rowMins(counts) > 0,]
dim(counts_filtered)
before_cd8_cxcl13_spots <- rownames(counts_filtered )
# add label
before@meta.data$cd8_cxcl13_label <- ifelse(rownames(before@meta.data) %in% before_cd8_cxcl13_spots,"CD8_CXCL13_Tcells","Others")
table(before@meta.data$cd8_cxcl13_label)

## Visualization
cols <- viridis(2,begin = 0.1,end = 0.6,option = "D",direction = -1)
names(cols) <- c("CD8_CXCL13_Tcells", "Others")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_0_cd8_cxcl13_Tcell_spatial_feature_plot.pdf",width = 12, height = 5)
SpatialDimPlot(after_ln,group.by = "cd8_cxcl13_label", pt.size.factor=1.3, alpha = c(1, 1), cols =cols)
SpatialDimPlot(after,group.by = "cd8_cxcl13_label", pt.size.factor=1.8, alpha = c(1, 1), cols =cols)
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_0_before_cd8_cxcl13_Tcell_spatial_feature_plot.pdf",width = 12, height = 5)
SpatialDimPlot(before,group.by = "cd8_cxcl13_label", pt.size.factor=2.3, alpha = c(1, 1), cols =cols)
dev.off()


########## 8. Check C1QC Macrophage location (FIg 4C,J) ######
after@active.assay <- "SCT"
after_ln@active.assay <- "SCT"
before@active.assay <- "SCT"
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/three_sample_C1QC_TAM_marker_spatial_featureplot.pdf")
SpatialFeaturePlot(before, features = c("CD68","LYZ",  "CD80","C1QC", "C1QB", "C1QA","CD44","CD74","MIF") , pt.size.factor=3, alpha = c(1, 1))
SpatialFeaturePlot(after,features = c("CD68","LYZ",  "CD80","C1QC", "C1QB", "C1QA","CD44","CD74","MIF") , pt.size.factor=1.8, alpha = c(1, 1))
SpatialFeaturePlot(after_ln, features = c("CD68","LYZ",  "CD80","C1QC", "C1QB", "C1QA","CD44","CD74","MIF") , pt.size.factor=1.3, alpha = c(1, 1))
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/three_sample_FCN1_TAM_marker_spatial_featureplot.pdf")
SpatialFeaturePlot(before, features = c("CD68","LYZ",  "FCN1","CD44","CD74","MIF") , pt.size.factor=3, alpha = c(1, 1))
SpatialFeaturePlot(after,features = c("CD68","LYZ",  "FCN1","CD44","CD74","MIF") , pt.size.factor=1.8, alpha = c(1, 1))
SpatialFeaturePlot(after_ln, features = c("CD68","LYZ",  "FCN1","CD44","CD74","MIF") , pt.size.factor=1.3, alpha = c(1, 1))
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/three_sample_CCL5_TAM_marker_spatial_featureplot.pdf")
SpatialFeaturePlot(before, features = c("CD68","LYZ",  "CCL5","CD44","CD74","MIF") , pt.size.factor=3, alpha = c(1, 1))
SpatialFeaturePlot(after,features = c("CD68","LYZ",  "CCL5","CD44","CD74","MIF") , pt.size.factor=1.8, alpha = c(1, 1))
SpatialFeaturePlot(after_ln, features = c("CD68","LYZ",  "CCL5","CD44","CD74","MIF") , pt.size.factor=1.3, alpha = c(1, 1))
dev.off()
## overlap LYZ C1QC CD68
## 1) after_ln
counts <- t(as.data.frame(after_ln@assays$Spatial@counts))
counts <- counts[,c("LYZ","CD68","C1QC")]
dim(counts)
# 3637    3
# All > 1
counts_filtered <- counts[rowMins(counts) > 0,]
#counts_filtered <- counts_filtered[rowSums(counts_filtered) >= 3,]
dim(counts_filtered)
# 1047   3
after_ln_C1QC_spots <- rownames(counts_filtered )
# add label
after_ln@meta.data$TAM_label <- ifelse(rownames(after_ln@meta.data) %in% after_ln_C1QC_spots,"C1QC+ TAM","Others")
table(after_ln@meta.data$TAM_label)

## 2) after
counts <- t(as.data.frame(after@assays$Spatial@counts))
counts <- counts[,c("LYZ","CD68","C1QC")]
dim(counts)
# 1244    3
# All > 1
counts_filtered <- counts[rowMins(counts) > 0,]
#counts_filtered <- counts_filtered[rowSums(counts_filtered) >= 3,]
dim(counts_filtered)
#  173   3
after_C1QC_spots <- rownames(counts_filtered )
# add label
after@meta.data$TAM_label <- ifelse(rownames(after@meta.data) %in% after_C1QC_spots,"C1QC+ TAM","Others")
table(after@meta.data$TAM_label)

## 3) before
counts <- t(as.data.frame(before@assays$Spatial@counts))
counts <- counts[,c("LYZ","CD68","C1QC")]
dim(counts)
# 1244    3
# All > 1
counts_filtered <- counts[rowMins(counts) > 1,]
#counts_filtered <- counts_filtered[rowSums(counts_filtered) >= 3,]
dim(counts_filtered)
#  173   3
before_C1QC_spots <- rownames(counts_filtered )
# add label
before@meta.data$TAM_label <- ifelse(rownames(before@meta.data) %in% before_C1QC_spots,"C1QC+ TAM","Others")
table(before@meta.data$TAM_label)

###### plot
cols <- viridis(2,begin = 0.1,end = 0.9,option = "E",direction = -1)
names(cols) <- c("C1QC+ TAM", "Others")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_0_C1QC_TAM_spatial_feature_plot.pdf",width = 12, height = 5)
SpatialDimPlot(after_ln,group.by = "TAM_label", pt.size.factor=1.3, alpha = c(1, 1), cols =cols)
SpatialDimPlot(after,group.by = "TAM_label", pt.size.factor=1.8, alpha = c(1, 1), cols =cols)
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/cutoff_1_C1QC_TAM_spatial_feature_plot.pdf",width = 12, height = 5)
SpatialDimPlot(before,group.by = "TAM_label", pt.size.factor=2.3, alpha = c(1, 1), cols =cols)
dev.off()

########## 9. Overlap ratio #######
### after
coexistence_celllabel <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after/matrix_spdata_plot_coexistence.rds")
coexistence_celllabel_ <- as.data.frame(coexistence_celllabel[,"co_existence"])
dim(coexistence_celllabel_)
rownames(coexistence_celllabel_) <- rownames(coexistence_celllabel)
labled_after <- AddMetaData(after,metadata = coexistence_celllabel_,col.name = "co_existence_label")
table(labled_after$co_existence_label)

after_overlap_ratio <- labled_after@meta.data %>% count(co_existence_label,label)


######### 10. CD8 CXCL13 vs ratio (Fig S11E,F,H,J)####
### after
cd8_counts <- after@assays$SCT@counts
cd8_counts <- as.data.frame(cd8_counts)
cd8_counts <- t(cd8_counts[c("CD8A","CXCL13","CD74","C1QC"),])
cd8_counts <- as.data.frame(cd8_counts)
#cd8_counts <- log10(cd8_counts)
cd8_counts$cell <- rownames(cd8_counts)

after@meta.data$cell <- rownames(after@meta.data)
after_merged_cd8_counts <- merge(cd8_counts, after@meta.data, by = "cell" )

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_no_log_ratio_cd8_boxtplot.pdf",width = 10)

ggboxplot(after_merged_cd8_counts,x = "label", y = "CD8A", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
        stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none") +
      
ggboxplot(after_merged_cd8_counts,x = "label", y = "CXCL13", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
      stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none")+

ggboxplot(after_merged_cd8_counts,x = "label", y = "CD74", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none")+

ggboxplot(after_merged_cd8_counts,x = "label", y = "C1QC", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none")

dev.off()
      
### after ln
cd8_counts <- after_ln@assays$SCT@counts
cd8_counts <- as.data.frame(cd8_counts)
cd8_counts <- t(cd8_counts[c("CD8A","CXCL13","CD74","C1QC"),])
cd8_counts <- as.data.frame(cd8_counts)
#cd8_counts <- log10(cd8_counts)
cd8_counts$cell <- rownames(cd8_counts)

after_ln@meta.data$cell <- rownames(after_ln@meta.data)
after_ln_merged_cd8_counts <- merge(cd8_counts, after_ln@meta.data, by = "cell" )

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_ln_no_log_ratio_cd8_boxtplot.pdf",width = 10)

ggboxplot(after_ln_merged_cd8_counts,x = "label", y = "CD8A", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) +theme(legend.position = "none") +
  
  ggboxplot(after_ln_merged_cd8_counts,x = "label", y = "CXCL13", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) +theme(legend.position = "none") +
  
  ggboxplot(after_ln_merged_cd8_counts,x = "label", y = "CD74", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none") +
  
  ggboxplot(after_ln_merged_cd8_counts,x = "label", y = "C1QC", color = "label") + scale_color_viridis(begin = 0.3, end = 0.9,discrete = T,option = "B") +
  stat_compare_means(method = "t.test",comparisons = list(c("low","medium"),c("medium","high"), c("low","high"))) + theme(legend.position = "none") 
dev.off()

######### 11. (After LN) SPATA for SiT  (Fig 5E)#########
#### after ln
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")

spata_obj <- transformSeuratToSpata(after_ln, sample_name = "after_ln",method="spatial", assay_name = "SCT")
spata_obj <- createTrajectories(spata_obj )

# get trajectory names 
getTrajectoryNames(object = spata_obj)
saveRDS(spata_obj,"~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_spata_obj.rds")
spata_obj <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_spata_obj.rds")

# visulization
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_ln_ratio_label.pdf")
plotTrajectory(object = spata_obj, 
               trajectory_name = "cd74_ratio",
               color_by = "label",
               #pt_clrp = "npg",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               display_image = FALSE) +  legendTop() 

plotTrajectory(object = spata_obj, 
               trajectory_name = "cd74_ratio",
               color_by = "seurat_clusters",
               smooth = TRUE, 
               pt_alpha = 0.25, 
               display_image = FALSE) +
  legendTop()

dev.off()

#  Trajectory trends via lineplots

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/SPATA_plotTrajectoryFeatures.pdf")

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio",
                       features = "nCount_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio",
                       features = "nFeature_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio",
                       features = "nCount_ISOG", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 


plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio",
                       features = "nCount_RATIO", 
                       smooth_method = "loess", 
                       smooth_span = 1, 
                       smooth_se = TRUE)

plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "cd74_ratio_new",
                               discrete_feature = "celltype", 
                               clrp = "inferno",
                               display_trajectory_parts = FALSE) 

plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "cd74_ratio",
                               discrete_feature = "label", 
                               clrp = "inferno",
                               display_trajectory_parts = FALSE) 
dev.off()


# gene-set names
genes_of_interest <- c("CD8A", "CXCL13", "CD74", "CCL4","CCL5", "CST7","NKG7",
                       "GZMK","GZMB","GZMA","MIF","IL1A","IL2","IL6","IL8")

# plot lineplot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_ln_SPATA_plotTrajectoryGenes.pdf")

plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "cd74_ratio", 
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE, 
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows 
)
plotTrajectoryHeatmap(object = spata_obj,
                      trajectory_name = "cd74_ratio", 
                      variables = genes_of_interest,show_rownames = T)
dev.off()

### pathway 
test <- c('CD8A',"CXCL13")
plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "cd74_ratio",
  gene_sets = "HM_WNT_BETA_CATENIN_SIGNALING",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "cd74_ratio",
  gene_sets = "HM_HYPOXIA",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()


##### after


######### 12. (After) SPATA for SiT  (Fig 5D) #########

spata_obj <- transformSeuratToSpata(after, sample_name = "after",method="spatial", assay_name = "SCT")
spata_obj <- createTrajectories(spata_obj )

# get trajectory names 
getTrajectoryNames(object = spata_obj)
saveRDS(spata_obj,"~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_spata_obj.rds")
spata_obj <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_spata_obj.rds")

# visulization
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_ln_SPATA_plotTrajectoryGenes_1.pdf")

plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "cd74_ratio", 
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE, 
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows 
)
plotTrajectoryHeatmap(object = spata_obj,
                      trajectory_name = "cd74_ratio", 
                      variables = genes_of_interest,show_rownames = T)
dev.off()
#  Trajectory trends via lineplots

pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_SPATA_plotTrajectoryFeatures_16.pdf")

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio_16",
                       features = "nCount_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio_16",
                       features = "nFeature_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio_16",
                       features = "nCount_ISOG", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 


plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "cd74_ratio_16",
                       features = "nCount_RATIO", 
                       smooth_method = "loess", 
                       smooth_span = 1, 
                       smooth_se = TRUE)

plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "cd74_ratio_16",
                               discrete_feature = "celltype", 
                               clrp = "inferno",
                               display_trajectory_parts = FALSE) 

plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "cd74_ratio_16",
                               discrete_feature = "label", 
                               clrp = "inferno",
                               display_trajectory_parts = FALSE) 
dev.off()


# gene-set names
genes_of_interest <- c("CD8A", "CXCL13", "CD74","CCL4", "CCL5", "GZMK")

# plot lineplot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/cd8_cxcl13/after_SPATA_plotTrajectoryGenes_16.pdf")

plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "cd74_ratio_16", 
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE, 
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows 
)
plotTrajectoryHeatmap(object = spata_obj,
                      trajectory_name = "cd74_ratio_16", 
                      variables = genes_of_interest,show_rownames = T)
dev.off()

### pathway 
test <- c('CD8A',"CXCL13")
plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "cd74_ratio",
  gene_sets = "HM_WNT_BETA_CATENIN_SIGNALING",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "cd74_ratio",
  gene_sets = "HM_HYPOXIA",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()


##### after



