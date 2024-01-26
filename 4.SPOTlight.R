######### 0. Libraries ##########
library(dplyr)
library(Seurat)
library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(NMF)
library(ggsci)
######### 1. Load data #########
### single cell data
scRNA_ori <-  readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")### spatial data

before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/3.Seurat/P1_before/P1_before.seurat.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after/p1_after.seurat.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after_LN/p1_after_LN.seurat.rds")

## change pre annotation data

#before <- readRDS("/thinker/tb269store/liuLab/Yin_yin/SiT/rawdata/p1_before.seurat.rds")
#after <- readRDS("/thinker/tb269store/liuLab/Yin_yin/SiT/rawdata/before_spdata.seurat.rds")
spdata <- after

######### 2. Preprocessing ########
#### 2.2 After
## delete FRC and endothelial cells
scRNA_ori <- subset(scRNA_ori, subset= orig.ident %in% c("Bcell", "Epithelia", "Fibroblast",
                                                         "Myeloid", "Pericytes", "Tcell") )
## change to single cell object
scRNA <- as.SingleCellExperiment(scRNA_ori)

### Feature selection
scRNA_feature_selected <- scuttle::logNormCounts(scRNA)

### Variance modelling
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(scRNA_feature_selected))
dec <- modelGeneVar(scRNA_feature_selected, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)
saveRDS(hvg,"~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_ldx_hvg.rds")
hvg <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_ldx_hvg.rds")


colLabels(scRNA_feature_selected) <- colData(scRNA_feature_selected)$orig.ident
# Compute marker genes
mgs <- scoreMarkers(scRNA_feature_selected, subset.row = genes)
# Keep relevent only
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
saveRDS(mgs_df,"~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_ldx_mgs_df.rds")
mgs_df <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/combined_sc_ldx_mgs_df.rds")
## no FRC and endothelial
mgs_df <- mgs_df[!mgs_df$cluster %in% c("FRC","Endothelial"), ]

## cell downsampling
idx <- split(seq(ncol(scRNA_feature_selected)), scRNA_feature_selected$label)
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
saveRDS(scRNA_feature_selected_ds, "~/onedrive/Work/phD/phd_project/SiT/results/sc/ldx_annotated_downsample_scRNA.rds")
scRNA_feature_selected_ds <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/ldx_annotated_downsample_scRNA.rds")

## no FRC and endothelial
other_cells <- rownames(scRNA_feature_selected_ds@colData)[!scRNA_feature_selected_ds@colData$label %in% c("FRC","Endothelial")]
scRNA_feature_selected_ds <- scRNA_feature_selected_ds[,other_cells]
                                  
######### SPOTlight Decomposition ##########
set.seed(1234)

res <- SPOTlight(
  x = scRNA_feature_selected_ds,
  y = as.matrix(before@assays$Spatial@counts),
  groups = as.character(scRNA_feature_selected_ds$label),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")


saveRDS(res, "~/onedrive/Work/phD/phd_project/SiT/results/before/before_spotlight_SiT_LDX.rds")
#res <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/spotlight_SiT_cd45_combined.rds")

# Extract deconvolution matrix
head(matrix_spdata <- res$mat)[, seq_len(3)]
mod <- res$NMF

######## Visualization #######
#setwd("/thinker/tb269store/liuLab/Yin_yin/SiT/")
setwd("~/onedrive/Work/phD/phd_project/SiT/")
#### 1. Topic profile
pdf("results/plotTopicProfiles.pdf")
plotTopicProfiles(
  x = mod,
  y = scRNA$annotation,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = scRNA$annotation,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)
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
#matrix_spdata[matrix_spdata < 0.] <- 0

# Define color palette
matrix_spdata_plot <- matrix_spdata
spdata <- before
#matrix_spdata_plot <- matrix_spdata[,colnames(matrix_spdata_plot) != "FRC" ]
pdf("results/before/downsample_100_plotSpatialScatterpie_no_FRC_Endo_all_1_2_3.pdf")
plotSpatialScatterpie(
  x = spdata,
  y = matrix_spdata_plot,
  cell_types = colnames(matrix_spdata_plot),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 1.2) +
  scale_fill_simpsons()
matrix_spdata_plot[matrix_spdata_plot < 0.1] <- 0
plotSpatialScatterpie(
  x = spdata,
  y = matrix_spdata_plot,
  cell_types = colnames(matrix_spdata_plot),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 1.2) +
  scale_fill_simpsons()

matrix_spdata_plot[matrix_spdata_plot < 0.2] <- 0
plotSpatialScatterpie(
  x = spdata,
  y = matrix_spdata_plot,
  cell_types = colnames(matrix_spdata_plot),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 1.2) +
  scale_fill_simpsons()
matrix_spdata_plot[matrix_spdata_plot < 0.3] <- 0
plotSpatialScatterpie(
  x = spdata,
  y = matrix_spdata_plot,
  cell_types = colnames(matrix_spdata_plot),
  #img = T,
  #slice = "image_png",
  scatterpie_alpha = 1,
  pie_scale = 1.2) +
  scale_fill_simpsons()
dev.off()
#### 4. Residuals (all zero)

spdata_ <- as.SingleCellExperiment(spdata)
spdata_$res_ss <- res[[2]][colnames(spdata_)]
xy <- GetTissueCoordinates(spdata)
spdata_$x <- xy[, 2]
spdata_$y <- -xy[, 1]
pdf("results/plotResiduals.pdf")
ggcells(spdata_, aes(x, y, color = res_ss)) +
  geom_point() +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_bw()
dev.off()
#### 5. 山峰图
## only look at one cell type
coor = spdata@images$image@coordinates[,c(5,4)]
st_meta <- cbind(coor,res$mat)
colnames(st_meta)[1:2] = c('x','y')
st_meta$x <- (st_meta$x - min(st_meta$x))/(max(st_meta$x) -
                                             min(st_meta$x))
st_meta$y <- (st_meta$y - min(st_meta$y))/(max(st_meta$y) -
                                             min(st_meta$y))
cellname <- colnames(st_meta)[-c(1:2)]
celltype <- "Tcell"

col_manual <- ggpubr::get_palette(palette = "lancet", 
                                  k = length(cellname))
col_manual <- pal_ucscgb()(length(cellname))
point_color <- col_manual[which(cellname == celltype)]

st_meta <- st_meta[,c("x","y", celltype)]

p <- ggplot(data = st_meta, ggplot2::aes(x, y)) + 
  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())

p <- p + ggplot2::geom_point(size = 1, color = point_color)

p <- p + ggplot2::stat_density2d(ggplot2::aes(color = ..level..), 
                                 size = 1)
p + scale_color_gradient(low = col_manual[1], high = col_manual[2])
#p + scale_color_gsea()
pdf("results/Tcell_plot.pdf")
print(p)
dev.off()

######## get results for stlearn ######
tmp = as.data.frame(res$mat)
row.names(tmp) <- row.names(spdata@meta.data)
write.csv(t(tmp[1:(length(tmp))]),"results/after_LN_deconvolution_result.csv")

###### check cluster percentages ########
## get deconvolution results
res <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before/before_spotlight_SiT_cd45_combined.rds")
tmp = as.data.frame(res$mat)
row.names(tmp) <- row.names(spdata@meta.data)
## combine into metadata
new_metadata <- cbind(spdata@meta.data,tmp)
cluster_cell_percentage <- new_metadata %>% group_by(seurat_clusters) %>% 
  summarise(avg_Bcell = median(Bcell),avg_Endothelial = median(Endothelial),
            avg_Epithelial = median(Epithelia),avg_Fibroblast = median(Fibroblast),
            avg_FRC = median(FRC),avg_Myeloid = median(Myeloid),
            avg_Pericytes = median(Pericytes),avg_Tcell = median(Tcell))
cluster_cell_percentage <- as.data.frame(cluster_cell_percentage)
colnames(cluster_cell_percentage)[2:9] <- rep("perc",8)
ggplot_cluster_cell_percentage <- rbind(cluster_cell_percentage[,c(1,2)], 
                                        cluster_cell_percentage[,c(1,3)],
                                        cluster_cell_percentage[,c(1,4)],
                                        cluster_cell_percentage[,c(1,5)],
                                        cluster_cell_percentage[,c(1,6)],
                                        cluster_cell_percentage[,c(1,7)],
                                        cluster_cell_percentage[,c(1,8)],
                                        cluster_cell_percentage[,c(1,9)])

ggplot_cluster_cell_percentage$celltype <- rep(colnames(tmp),1,each =6)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/after_LN/celltype_assigning_6_cluster.pdf")
ggplot(data = ggplot_cluster_cell_percentage, aes(x = seurat_clusters, y = perc,fill = celltype)) + 
  geom_bar( stat = "identity", position = "stack")
dev.off()
##cluster 分得不好
