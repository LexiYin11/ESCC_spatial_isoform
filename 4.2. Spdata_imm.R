library(SeuratData)
library(SpatialExperiment)
library(dplyr)
library(Seurat)
library(patchwork)
####### 0. Data import ########
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/3.Seurat/P1_before/P1_before.seurat.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after/p1_after.seurat.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after_LN/p1_after_LN.seurat.rds")

scRNA <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/sce_anno_res_0.2.rds")

####### 1.0 cell annotation #######
### after LN 
after_ln@meta.data$celltype <- ifelse(after_ln$seurat_clusters == "0","Bcells, Myeloid and Tcells RGS1+",
                                      ifelse(after_ln$seurat_clusters == "1","Epithelial cells COL1A2+",
                                             ifelse(after_ln$seurat_clusters == "2","Bcells and Fibroblast IGLC2+",
                                                    ifelse(after_ln$seurat_clusters == "3","Bcells , Myeloid and Tcell APOC1+",
                                                           ifelse(after_ln$seurat_clusters == "4","Epithelial cells CXCL13+",
                                                                  ifelse(after_ln$seurat_clusters == "5","Epithelial cells FABP5+","oh no"))))))
table(after_ln$seurat_clusters)
table(after_ln$celltype)
### after 
after@meta.data$celltype <- ifelse(after$seurat_clusters == "0",'FRC and Epithelial COL1A2+',
                                   ifelse(after$seurat_clusters == "1",'Fibroblast, Pericytes and Endothelial MYH11+',
                                          ifelse(after$seurat_clusters == "2",'Bcell, Tcell, Myeloid IGHG4+',
                                                 ifelse(after$seurat_clusters == "3",'Fibroblast, Pericytes and Endothelial COL1A1-',
                                                        ifelse(after$seurat_clusters == "4",'Epithelial CRISP3+',
                                                               ifelse(after$seurat_clusters == "5",'FRC and Epithelial FABP5+',
                                                                      ifelse(after$seurat_clusters == "6",'Bcell, Tcell, Myeloid CXCL13+',
                                                                             ifelse(after$seurat_clusters == "7",'Epithelial MUC4+',"no"))))))))


table(after$seurat_clusters)
table(after$celltype)

## before
before@meta.data$celltype <- ifelse(before$seurat_clusters == "0",'FRC and Epithelial KRT178-',
                                    ifelse(before$seurat_clusters == "1",'FRC,Fibroblast and Epithelial CXCL14+',
                                           ifelse(before$seurat_clusters == "2",'FRC and Epithelial KRT4+',
                                                  ifelse(before$seurat_clusters == "3",'Bcell, FRC and Epithelial IFITM3+',
                                                         ifelse(before$seurat_clusters == "4",'FRC and Epithelial IGFBP5+',
                                                                ifelse(before$seurat_clusters == "5",'Epithelial MUC4+',"no"))))))


######  2.1 extract imm cell only ######
#### before ####
before_imm_cells <- subset(before, subset= celltype %in% c('Epithelial cells (CSTB+), T cells, B cells, Macrophage',
                                                          'Bcells T cells',
                                                         'Tcell Macrophage Fibroblast '))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/before/only_imm_spatialDimPlot.pdf")
p1 <- SpatialDimPlot(before_imm_cells,group.by = "celltype", pt.size.factor = 6,label = F, label.size = 4)
p2 <- SpatialDimPlot(before_imm_cells, crop= F,pt.size.factor = 1,label = F, label.size = 24)
p1
p2
dev.off()
#### after ####
after_imm_cells <- subset(after, subset= celltype %in% c('Malignant Epithelial cells , Macrophage, T cell','B cells, Tcells, Fibroblast,Macrophage',
                                                         'Tcell Macrophage Fibroblast '))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/after/only_imm_spatialDimPlot.pdf")
p1 <- SpatialDimPlot(after_imm_cells,group.by = "celltype", pt.size.factor = 2.5,label = F, label.size = 4)
p2 <- SpatialDimPlot(after_imm_cells, crop= F,pt.size.factor = 1,label = F, label.size = 2)
p1 + p2
dev.off()
#### after_LN ####
after_LN_imm_cells <- after_ln
pdf("~/onedrive/Work/phD/phd_project/SiT/results/after_LN/only_imm_spatialDimPlot.pdf")
p1 <- SpatialDimPlot(after_LN_imm_cells,group.by = "celltype", pt.size.factor = 1.5,label = TRUE, label.size = 4)
p2 <- SpatialDimPlot(after_LN_imm_cells, crop= F,pt.size.factor = 1,label = TRUE, label.size = 2)
p1 + p2
dev.off()

######  2.2 remove cell marker of all cell types except T cells #####
markers <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/res_0.2_cosg_cluster_marker.csv")
only_T_cells_markers <- markers[,c(1,2,6,7)]
only_T_cells_markers <- c(unlist(c(only_T_cells_markers)))
## before
before_imm_cells_T_genes <- before_imm_cells
before_imm_cells_T_genes@assays$Spatial@counts <- before_imm_cells_T_genes@assays$Spatial@counts[rownames(before_imm_cells_T_genes@assays$Spatial@counts) %in% only_T_cells_markers, ]

## after
after_imm_cells_T_genes <- after_imm_cells
after_imm_cells_T_genes@assays$Spatial@counts <- after_imm_cells_T_genes@assays$Spatial@counts[rownames(after_imm_cells_T_genes@assays$Spatial@counts) %in% only_T_cells_markers, ]

## after
after_LN_imm_cells_T_genes <- after_LN_imm_cells
after_LN_imm_cells_T_genes@assays$Spatial@counts <- after_LN_imm_cells_T_genes@assays$Spatial@counts[rownames(after_LN_imm_cells_T_genes@assays$Spatial@counts) %in% only_T_cells_markers, ]

####### 2.3  reduction & recluster  ######
## before
before_imm_cells_T_genes_normalized <- SCTransform(before_imm_cells_T_genes, assay = "Spatial", verbose = FALSE)
before_imm_cells_T_genes_normalized <- RunPCA(before_imm_cells_T_genes_normalized, assay = "SCT", verbose = FALSE)
before_imm_cells_T_genes_normalized <- FindNeighbors(before_imm_cells_T_genes_normalized, reduction = "pca", dims = 1:30)
before_imm_cells_T_genes_clustered <- FindClusters(before_imm_cells_T_genes_normalized, verbose = FALSE)
before_imm_cells_T_genes_clustered_umap <- RunUMAP(before_imm_cells_T_genes_clustered, reduction = "pca", dims = 1:30)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/before/imm_cells_T_genes_SpatialDimPlot.pdf")
SpatialDimPlot(before_imm_cells_T_genes_clustered_umap, crop= F,pt.size.factor = 1,label = TRUE, label.size = 4)
dev.off()
table(before_imm_cells_T_genes_clustered_umap$seurat_clusters)

## after_LN
after_LN_imm_cells_T_genes_normalized <- SCTransform(after_LN_imm_cells_T_genes, assay = "Spatial", verbose = FALSE)
after_LN_imm_cells_T_genes_normalized <- RunPCA(after_LN_imm_cells_T_genes_normalized, assay = "SCT", verbose = FALSE)
after_LN_imm_cells_T_genes_normalized <- FindNeighbors(after_LN_imm_cells_T_genes_normalized, reduction = "pca", dims = 1:30)
after_LN_imm_cells_T_genes_clustered <- FindClusters(after_LN_imm_cells_T_genes_normalized, verbose = FALSE)
after_LN_imm_cells_T_genes_clustered_umap <- RunUMAP(after_LN_imm_cells_T_genes_clustered, reduction = "pca", dims = 1:30)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/after_LN/imm_cells_T_genes_SpatialDimPlot.pdf")
SpatialDimPlot(after_LN_imm_cells_T_genes_clustered_umap, crop= F,pt.size.factor = 1,label = TRUE, label.size = 2)
dev.off()
table(after_LN_imm_cells_T_genes_clustered_umap$seurat_clusters)

## after
after_imm_cells_T_genes_normalized <- SCTransform(after_imm_cells_T_genes, assay = "Spatial", verbose = FALSE)
after_imm_cells_T_genes_normalized <- RunPCA(after_imm_cells_T_genes_normalized, assay = "SCT", verbose = FALSE)
after_imm_cells_T_genes_normalized <- FindNeighbors(after_imm_cells_T_genes_normalized, reduction = "pca", dims = 1:30)
after_imm_cells_T_genes_clustered <- FindClusters(after_imm_cells_T_genes_normalized, verbose = FALSE)
after_imm_cells_T_genes_clustered_umap <- RunUMAP(after_imm_cells_T_genes_clustered, reduction = "pca", dims = 1:30)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/after/imm_cells_T_genes_SpatialDimPlot.pdf")
SpatialDimPlot(after_imm_cells_T_genes_clustered_umap, crop= F,pt.size.factor = 1.3,label = TRUE, label.size = 2)
dev.off()
table(after_imm_cells_T_genes_clustered_umap$seurat_clusters)

### extract only 0,2,3
#after_imm_cells<- subset(after_imm_cells_clustered_umap, subset= seurat_clusters %in% c(0,2,3))

######## 2.4 extract markers #####
### before
before_imm_cells_T_genes_allmarkers <- FindAllMarkers(before_imm_cells_T_genes_clustered_umap)
before_imm_cells_T_genes_allmarkers.markers = before_imm_cells_T_genes_allmarkers %>% select(gene, everything()) %>% subset(p_val<0.05)
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
before_imm_cells_T_genes_allmarkers_top100 = before_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
before_imm_cells_T_genes_allmarkers_top10 = before_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(before_imm_cells_T_genes_allmarkers_top100, "~/onedrive/Work/phD/phd_project/SiT/results/before/findallmarker_cluster_marker_top100.csv", row.names = F)
write.csv(before_imm_cells_T_genes_allmarkers_top10, "~/onedrive/Work/phD/phd_project/SiT/results/before/findallmarker_cluster_marker_top10.csv", row.names = F)

### after
after_imm_cells_T_genes_allmarkers <- FindAllMarkers(after_imm_cells_T_genes_clustered_umap)
after_imm_cells_T_genes_allmarkers.markers = after_imm_cells_T_genes_allmarkers %>% select(gene, everything()) %>% subset(p_val<0.05)
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
after_imm_cells_T_genes_allmarkers_top100 = after_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
after_imm_cells_T_genes_allmarkers_top10 = after_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(after_imm_cells_T_genes_allmarkers_top100, "~/onedrive/Work/phD/phd_project/SiT/results/after/findallmarker_cluster_marker_top100.csv", row.names = F)
write.csv(after_imm_cells_T_genes_allmarkers_top10, "~/onedrive/Work/phD/phd_project/SiT/results/after/findallmarker_cluster_marker_top10.csv", row.names = F)

### after_LN
after_LN_imm_cells_T_genes_allmarkers <- FindAllMarkers(after_LN_imm_cells_T_genes_clustered_umap)
after_LN_imm_cells_T_genes_allmarkers.markers = after_LN_imm_cells_T_genes_allmarkers %>% select(gene, everything()) %>% subset(p_val<0.05)
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
after_LN_imm_cells_T_genes_allmarkers_top100 = after_LN_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
after_LN_imm_cells_T_genes_allmarkers_top10 = after_LN_imm_cells_T_genes_allmarkers.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(after_LN_imm_cells_T_genes_allmarkers_top100, "~/onedrive/Work/phD/phd_project/SiT/results/after_LN/findallmarker_cluster_marker_top100.csv", row.names = F)
write.csv(after_LN_imm_cells_T_genes_allmarkers_top10, "~/onedrive/Work/phD/phd_project/SiT/results/after_LN/findallmarker_cluster_marker_top10.csv", row.names = F)

######## 3.1 scRNA re-annotation Tcells only ####
anno_table <- data.frame("cluster" = c(0,1:10),"celltype" = c("CD8+_T_cells_effector",
                                                              "CD4+_T_cells", "B_cells_plasma",
                                                              "Macrophage",
                                                              "B_cells_follicular","CD8+_T_cells_terminal_effector",
                                                              "CD8+_T cells_Naive","Basophils",
                                                              "Endothelial_cells",
                                                              "B_cells_plasma", "Smooth_muscle_cells"))

scRNA@meta.data$barcode =rownames(scRNA@meta.data)
scRNA@meta.data=merge(scRNA@meta.data,anno_table,by.x="seurat_clusters",by.y="cluster")
rownames(scRNA@meta.data)= scRNA@meta.data$barcode
table(scRNA@meta.data$celltype.y)
table(scRNA@meta.data$seurat_clusters)
## only keep T cells
scRNA_Tcell <- subset(scRNA,subset = celltype.y %in% c("CD4+_T_cells","CD8+_T_cells_effector","CD8+_T cells_Naive","CD8+_T_cells_terminal_effector"))
table(scRNA_Tcell@meta.data$celltype.y)   
scRNA <- scRNA_Tcell
## change to single cell object
scRNA <- as.SingleCellExperiment(scRNA)
markers <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/res_0.2_cosg_cluster_marker.csv")
markers <- mar
### only T cells
only_T_cells_markers
markers_cluster <- data.frame("cluster" = rep(c(0,1,5,6),each = 100), markers = unlist(c(only_T_cells_markers)))
cluster_celltype <- unique(scRNA@meta.data[,c("seurat_clusters","celltype.y")])
merged_markers_cluster <- merge(markers_cluster, cluster_celltype,by.x ="cluster", by.y  = "seurat_clusters")
scRNA_feature_selected <- scuttle::logNormCounts(scRNA)

### Variance modelling
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(scRNA_feature_selected))
dec <- modelGeneVar(scRNA_feature_selected, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)
colLabels(scRNA_feature_selected) <- colData(scRNA_feature_selected)$celltype.y
# Compute marker genes
mgs <- scoreMarkers(scRNA_feature_selected, subset.row = genes)
# Keep relevent only
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.6, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)


## cell downsampling
idx <- split(seq(ncol(scRNA_feature_selected)), scRNA_feature_selected$celltype.y)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
scRNA_feature_selected_ds <- scRNA_feature_selected[, unlist(cs_keep)]
setwd("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/")
#saveRDS(scRNA_feature_selected_ds, "annotated_downsample_scRNA.rds")
#scRNA_ds <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/CD45_combined_hpca_main_annotated_downsample_scRNA.rds")
## 结果一般
######## 3.2 SPOTlight Decomposition ##########
set.seed(1234)
#scRNA <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_new_downsampling_sce.rds")## add celltypes
#mgs_df <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_new_mgs_df.rds")
#hvg <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_new_hvg.rds")
#mgs_df <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_new_mgs_df.rds")
#hvg <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_new_hvg.rds")
res <- SPOTlight(
  x = scRNA_feature_selected_ds,
  y = as.matrix(after_LN_imm_cells_T_genes_clustered_umap@assays$SCT@counts),
  groups = as.character(scRNA_feature_selected_ds$celltype.y),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

saveRDS(res, "~/onedrive/Work/phD/phd_project/SiT/results/after/after_Tcells_spotlight_SiT_self.rds")
#res <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/spotlight_SiT_cd45_combined.rds")
saveRDS(res, "~/onedrive/Work/phD/phd_project/SiT/results/before/before_Tcells_spotlight_SiT_self.rds")
saveRDS(res, "~/onedrive/Work/phD/phd_project/SiT/results/after_LN/after_LN_Tcells_spotlight_SiT_self.rds")

# Extract deconvolution matrix
head(matrix_spdata <- res$mat)[, seq_len(3)]
mod <- res$NMF

######## 3.3 Deconvlution Visualization #######
#setwd("/thinker/tb269store/liuLab/Yin_yin/SiT/")
setwd("~/onedrive/Work/phD/phd_project/SiT/")
#### 1. Topic profile
pdf("results/after/plotTopicProfiles.pdf")
plotTopicProfiles(
  x = mod,
  y = scRNA_feature_selected_ds$celltype.y,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

dev.off()
#### 2. Correlation Matrix
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
pdf("results/plotInteractions.pdf")
plotCorrelationMatrix(matrix_spdata)
plotInteractions(matrix_spdata, which = "heatmap", metric = "prop")
plotInteractions(matrix_spdata, which = "heatmap", metric = "jaccard")
plotInteractions(matrix_spdata, which = "network")
dev.off()
#### 3. Scatterpie
## you can choose interest gene you want show
ct <- colnames(matrix_spdata)
matrix_spdata_pie <- matrix_spdata
# Define color palette
pdf("results/after_LN/Tcells_plotSpatialScatterpie.pdf")
plotSpatialScatterpie(
  x = after_LN_imm_cells_T_genes_clustered_umap,
  y = matrix_spdata_pie,
  cell_types = colnames(matrix_spdata_pie),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_jco()
matrix_spdata_pie[matrix_spdata_pie < 0.2] <- 0
plotSpatialScatterpie(
  x = after_LN_imm_cells_T_genes_clustered_umap,
  y = matrix_spdata_pie,
  cell_types = colnames(matrix_spdata_pie),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_jco()

matrix_spdata_pie[matrix_spdata_pie < 0.4] <- 0
plotSpatialScatterpie(
  x = after_LN_imm_cells_T_genes_clustered_umap,
  y = matrix_spdata_pie,
  cell_types = colnames(matrix_spdata_pie),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_jco()
matrix_spdata_pie[matrix_spdata_pie < 0.5] <- 0
plotSpatialScatterpie(
  x = after_LN_imm_cells_T_genes_clustered_umap,
  y = matrix_spdata_pie,
  cell_types = colnames(matrix_spdata_pie),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_jco()
dev.off()


#### 4. Check feature #######
##### before ####
before_imm_cells_T_genes_clustered_umap@active.assay <- "SCT"
marker_list <- c("CD4","CD8A","CD3D","CD3E","CD27","CCR7","GZMK","IL7R","CD44","TIGIT")

p1 = SpatialFeaturePlot(before_imm_cells,features = "CD3D", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p2 = SpatialFeaturePlot(before_imm_cells,features = "CD8A", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p3 = SpatialFeaturePlot(before_imm_cells,features = "CD4", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p4 = SpatialFeaturePlot(before_imm_cells,features = "CD27", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p5 = SpatialFeaturePlot(before_imm_cells,features = "CCR7", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p6 = SpatialFeaturePlot(before_imm_cells,features = "GZMK", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p7 = SpatialFeaturePlot(before_imm_cells,features = "IL7R", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p8 = SpatialFeaturePlot(before_imm_cells,features = "CD44", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p9 = SpatialFeaturePlot(before_imm_cells,features = "TIGIT", pt.size.factor = 5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/before/only_T_cell_markers_SpatialFeaturePlot.pdf",width = 10,height = 15)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)
dev.off()

##### after ####
after_imm_cells_T_genes_clustered_umap@active.assay <- "SCT"
marker_list <- c("CD4","CD8A","CD3D","CD3E","CD27","CCR7","GZMK","IL7R","CD44","TIGIT")

p1 = SpatialFeaturePlot(after_imm_cells,features = "CD3D", pt.size.factor = 2.5) +
scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p2 = SpatialFeaturePlot(after_imm_cells,features = "CD8A", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p3 = SpatialFeaturePlot(after_imm_cells,features = "CD4", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p4 = SpatialFeaturePlot(after_imm_cells,features = "CD27", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p5 = SpatialFeaturePlot(after_imm_cells,features = "CCR7", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p6 = SpatialFeaturePlot(after_imm_cells,features = "GZMK", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p7 = SpatialFeaturePlot(after_imm_cells,features = "IL7R", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p8 = SpatialFeaturePlot(after_imm_cells,features = "CD44", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
p9 = SpatialFeaturePlot(after_imm_cells,features = "TIGIT", pt.size.factor = 2.5) +
  scale_fill_gradientn(limits = c(0,3),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/after/only_T_cell_markers_SpatialFeaturePlot.pdf",width = 10,height = 15)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)
dev.off()



######### 1. Heatmap #####
## top marker vs cluster
### #1 plot
c(unlist(c(marker_cosg$names[1,])))
marker_cosg <- read.table("/public/home/ac8q5y82z9/data/3.Results/cell_identity/marker_cosg_qc4_set2_harmony_old.csv",header = T)
DoHeatmap(scRNA,features = c(unlist(c(marker_cosg$names[1:5,]))),label = F,slot = "scale.data")
ggsave("/public/home/ac8q5y82z9/data/3.Results/cell_identity/top5_markers,heatmap.pdf",device = "pdf",width = 50,height = 40,units = "cm")

## #2 plot
library(pheatmap)
colanno=mye.seu@meta.data[,c("CB","celltype")]
colanno=colanno%>%arrange(celltype)
rownames(colanno)=colanno$CB
colanno$CB=NULL
colanno$celltype=factor(colanno$celltype,levels = unique(colanno$celltype))

rowanno=markerdf1
rowanno=rowanno%>%arrange(celltype)

mat4=mye.seu[["RNA"]]@scale.data[rowanno$gene,rownames(colanno)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！

pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno,
         gaps_row=as.numeric(cumsum(table(rowanno$celltype))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]),
         filename="heatmap.2.pdf",width=11,height = 7
)

######### 2. Dotplot ########
### top marker (few) expression level and percent expressed vs cluster
marker <- c("ACKR1","RAMP2",
            "LUM","COL3A1",
            "KRT14","KRT5",
            "CD69","CD52")
DotPlot(scedata, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

####### 3.Featureplot ######
## marker location in all clusers vs expression
T.cell.marker = c("CD3D",'CD3E','CD2')

Fib.cell.marker = c('COL1A1','DCN','C1R')
Myeioid.cell.marker = c('LYZ','CD68','TYROBP')
B.cell.marker= c('CD79A','MZB1','MS4A1')
Endothelial.cell.marker= c('CLDN5','FLT1','RAMP2')
Mast.cell.marker= c('CPA3','TPSAB1','TPSB2')
DC.cell.marker= c('LILRA4','CXCR3','IRF7')
all_markers <- c(T.cell.marker,Fib.cell.marker,Myeioid.cell.marker,B.cell.marker,
                 Endothelial.cell.marker, Mast.cell.marker,DC.cell.marker   ) 
scRNA <- sce_harmony
#imm_markers <- c("CD2","CD3E","CD3D","CD8A","FOXP3","NKG7","CD19","JCHAIN")
pdf("/public/home/ac8q5y82z9/data/3.Results/Markers/Tcell_markers_feature_plot_qc6.pdf")
cols = c("gray", "coral2")
FeaturePlot(scRNA, features = T.cell.marker,cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/Markers/Myeioid_cell_markers_feature_plot_qc6.pdf")
cols = c("gray", "coral2")
FeaturePlot(scRNA, features = Myeioid.cell.marker,cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/Markers/B_cell_markers_feature_plot_qc6.pdf")
cols = c("gray", "coral2")
FeaturePlot(scRNA, features = B.cell.marker,cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/Markers/all_cell_markers_feature_plot_qc6.pdf")
cols = c("gray", "coral2")
FeaturePlot(scRNA, features = all_markers,cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/Markers/CD8_CD4_markers_feature_plot_qc6.pdf")
cols = c("gray", "coral2")
FeaturePlot(scRNA, features = c("CD4","CD8A","CD3B","CD3D","CD3E","CD3G","NCAM1","CD56","KLRB1","NKG7","PTPRC"
),cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()


