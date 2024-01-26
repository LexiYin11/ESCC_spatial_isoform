########### 0. Library import ########
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
RenameGenesSeurat <- function(obj = obj) {
  ## Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  #print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  newnames <- toupper(rownames(RNA))
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

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
########### 1. Create object ##########
setwd("/public/home/ac8q5y82z9/data/2.pre_Seurat/SiT")
folders=list.files('./')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
for (i in 1:3){
  sceList[[i]] <- RenameGenesSeurat(sceList[[i]])
}
sample_list <- c("after", "after_ln", "before")

## merged
sce_merged <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]]),
                    add.cell.ids = folders, 
                    project = "SiT_scRNA")
saveRDS(sce_merged,"/public/home/ac8q5y82z9/data/3.Results/SiT/merged/sce_merged.rds")
sce_merged <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/sce_merged.rds")
########### 2. QC - before ##############
sce_qc <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/sce_merged.rds")
sce_qc <- sce_merged
DefaultAssay(sce_qc) <- "RNA"
sce_qc@meta.data$project <- "SiT_merged"
##计算线粒体和红细胞基因比例
sce_qc[["percent.mt"]] <- PercentageFeatureSet(sce_qc, pattern = "^MT-")
#}  
##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
setwd("/public/home/ac8q5y82z9/data/3.Results/SiT/merged")
col.num <- length(levels(as.factor(sce_qc@meta.data$orig.ident)))
violin <-VlnPlot(sce_qc, group.by = "project",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 cols =rainbow(col.num),
                 pt.size = 0, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste0("QC_pre/merged_vlnplot_before_qc.pdf"), plot = violin, width = 12, height = 6) 
ggsave(paste0("QC_pre/merged_vlnplot_before_qc.png"), plot = violin, width = 12, height = 6) 

### check T cell
T.cell.marker = c("CD3D",'CD3E',"CD3G")
B.cell.marker= c('CD79A','MZB1','MS4A1')
Myeioid.cell.marker = c('LYZ','CD68','TYROBP')
TBandMyeloidcell <- c(T.cell.marker, B.cell.marker,Myeioid.cell.marker)

only_T_cell <- sce_qc[T.cell.marker,]
only_T_cell <- only_T_cell[,colSums(only_T_cell) > 0]

only_T_cell

## epithelial cell
epi.marker <- c("EPCAM","KRT5")

epi_cell <- sce_qc[epi.marker,]
epi_cell <- epi_cell[,colSums(epi_cell) > 0]

epi_cell
########### 3. QC - after ##############
#### qc 1 
minGene=500
maxGene=4000
pctMT=5

sce_qc1 <-sce_qc
sce_qc1@meta.data$project <- "SiT"
##计算线粒体和红细胞基因比例
sce_qc1[["percent.mt"]] <- PercentageFeatureSet(sce_qc1, pattern = "^MT-")
sce_qc1 <- subset(sce_qc1, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(sce_qc1@meta.data$orig.ident)))
#sceList[[i]] <- sce_qc1
violin <- VlnPlot(sce_qc1, group.by = "project",
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  cols =rainbow(col.num), 
                  pt.size = 0, 
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC1/after_vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC1/after_vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

## check T cell
only_T_cell <- sce_qc1[T.cell.marker,]
only_T_cell <- only_T_cell[,colSums(only_T_cell) > 0]

only_T_cell

########### 4. Normalized and scale and PCA ######
sce_qc <- sce_qc1
sce_qc_normalized <- NormalizeData(sce_qc, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce_qc_feature <- FindVariableFeatures(sce_qc_normalized)
sce_qc_scaled <- ScaleData(sce_qc_feature )
sce_qc_PCA <- RunPCA(sce_qc_scaled, features = VariableFeatures(object = sce_qc_scaled))

## harmony
sce_harmony <- RunHarmony(sce_qc_PCA , group.by.vars = "orig.ident") 

########### 5. PC number #####

# Determine percent of variation associated with each PC
pct <- sce_harmony [["pca"]]@stdev / sum( sce_harmony [["pca"]]@stdev) * 100


# Calculate cumulative percents for each PC
cumu <- cumsum(pct)


# 1. Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
## qc3 41
## qx4 41
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1


# 2. last point where change of % of variation is more than 0.1%.
co2
## qc4 16

# 3. Minimum of the two calculation
pcs <- min(co1, co2)
pcs

########### 6. UMAP and TSNE and cluster #########

sce_harmony_umap <- sce_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:17) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:17) %>%
  FindClusters( dims = 1:17,resolution = 0.2) %>% 
  identity()
sce_harmony_tsne <- sce_harmony_umap %>% 
  RunTSNE(reduction = "harmony", dims = 1:17)

saveRDS(sce_harmony_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/sce_harmony_old_clustered_qc1_pc17_res_0.5/res_0.5.rds")
sce_harmony_tsne <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/sce_harmony_old_clustered_qc1_pc17_res_0.5/res_0.5.rds")


########### 7. check doublelet #####
## pK Identification
sweep.res.list_sce_harmony <- paramSweep_v3(sce_harmony_tsne, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_sce_harmony)
sweep.stats_sce_harmony <- summarizeSweep(sweep.res.list_sce_harmony, GT = FALSE)
bcmvn_sce_harmony <- find.pK(sweep.stats_sce_harmony) 
## BEST CUTOFF：0.005
mpK<-as.numeric(as.vector(bcmvn_sce_harmony$pK[which.max(bcmvn_sce_harmony$BCmetric)]))

## (3) Homotypic Doublet Proportion Estimate 
annotations <- sce_harmony_tsne@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(sce_harmony_tsne)*8*1e-6 
DoubletRate = 0.224392
nExp_poi <- round(DoubletRate*length(sce_harmony_tsne$seurat_clusters)) 
# Ratio
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 5109

## (4) Run DoubletFinder with varying classification stringencies 
sce_harmony_tsne <- doubletFinder_v3(sce_harmony_tsne, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
sce_harmony_tsne <- doubletFinder_v3(sce_harmony_tsne, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)

### visualization
sce_harmony_tsne@meta.data[,"DF_hi.lo"] <- sce_harmony_tsne@meta.data$DF.classifications_0.25_0.005_6294
sce_harmony_tsne@meta.data$DF_hi.lo[which(sce_harmony_tsne@meta.data$DF_hi.lo == "Doublet" & sce_harmony_tsne@meta.data$DF.classifications_0.25_0.005_5109 == "Singlet")] <- "Doublet-Low Confidience"
sce_harmony_tsne@meta.data$DF_hi.lo[which(sce_harmony_tsne@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(sce_harmony_tsne@meta.data$DF_hi.lo)
#Doublet-High Confidience  Doublet-Low Confidience                  Singlet 
#5109                     1185                    21755

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/2_doubletFinder.pdf")
DimPlot(sce_harmony_tsne, reduction = "tsne", group.by ="DF_hi.lo",cols =c("black","red","gold"))
dev.off()

## remove double let
sce_remove_doublet <- subset(sce_harmony_tsne, subset = DF_hi.lo == "Singlet")

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/after_2_doubletFinder.pdf")
DimPlot(sce_remove_doublet, reduction = "tsne", group.by ="DF_hi.lo",cols =c("black","red","gold"))
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/removed_doublet_dimplot.pdf")
DimPlot(sce_remove_doublet, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.5, raster = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
dev.off()
saveRDS(sce_remove_doublet, file = "/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/sc_sce_harmony_after_remove_doublet.rds")
sce_remove_doublet <- readRDS("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/sc_sce_harmony_after_remove_doublet.rds")

########### 8 Reclustering ########

sce_remove_doublet_normalized <- NormalizeData(sce_remove_doublet, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce_remove_doublet_feature <- FindVariableFeatures(sce_remove_doublet_normalized)
sce_remove_doublet_scaled <- ScaleData(sce_remove_doublet_feature )
sce_remove_doublet_PCA <- RunPCA(sce_remove_doublet_scaled, features = VariableFeatures(object = sce_remove_doublet_scaled))

pct <- sce_remove_doublet[["pca"]]@stdev / sum( sce_remove_doublet[["pca"]]@stdev) * 100


# Calculate cumulative percents for each PC
cumu <- cumsum(pct)


# 1. Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
## qc3 41
## qx4 41
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1


# 2. last point where change of % of variation is more than 0.1%.
co2
## qc4 16

# 3. Minimum of the two calculation
pcs <- min(co1, co2)
pcs

sce_remove_doublet_umap <- sce_remove_doublet%>% 
  RunUMAP(reduction = "harmony", dims = 1:17) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:17) %>%
  FindClusters( dims = 1:17,resolution = 0.5) %>% 
  identity()
sce_remove_doublet_tsne <- sce_remove_doublet_umap %>% 
  RunTSNE(reduction = "harmony", dims = 1:17)

saveRDS(sce_remove_doublet_tsne,"/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_sce_remove_doublet_tsne_clustered.rds")

########### 9. Plot ########
p1 <- DimPlot(sce_remove_doublet_tsne, reduction = "tsne", group.by = "orig.ident", 
              pt.size=0.5, raster = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p2 <- DimPlot(sce_remove_doublet_tsne, reduction = "tsne", group.by = "seurat_clusters", 
              pt.size=0.5, raster = F, cols = celltype_cols)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/dimplot_harmony_qc1_pc_17.pdf")
p1
p2

dev.off()

########### 10. Markers  ######
sce_cluster <- sce_remove_doublet_tsne
marker_cosg <- cosg(
  sce_cluster,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)
write.csv(marker_cosg$names, paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_cosg_cluster_marker.csv"), row.names = F)
write.csv(marker_cosg$names[1:10,], paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_cosg_cluster_marker_top10.csv"), row.names = F)

all.markers <- FindAllMarkers(sce_cluster)
all.markers <- markers
all.markers = all.markers%>% select(gene, everything()) %>% subset(p_val<0.05)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/findallmarker_cluster_marker.csv", row.names = F)
write.csv(top10, paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/top10_findallmarkercluster_marker_top10.csv"), row.names = F)

########### 11. Annotation prep ######
#### HumanPrimaryCellAtlasData
#refdata1 <- HumanPrimaryCellAtlasData()
#saveRDS(refdata1,"~/onedrive/Work/phD/phd_project/Gender_imm/results/HumanPrimaryCellAtlasData.rds")
refdata1 <- readRDS("/public/home/ac8q5y82z9/data/3.Results/cell_identity/HumanPrimaryCellAtlasData.rds")

testdata1 <- GetAssayData(sce_cluster, slot="data")
clusters <- sce_cluster@meta.data$seurat_clusters

cellpred_main <- SingleR(test = testdata1, ref = refdata1, labels = refdata1$label.main, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype=cellpred_main$labels, stringsAsFactors = F)
write.csv(celltype_main,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_HumanPrimaryCellAtlasData_celltype_singleR.csv"),row.names = F)

cellpred_fine <- SingleR(test = testdata1, ref = refdata1, labels = refdata1$label.fine, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype=cellpred_fine$labels, stringsAsFactors = F)
write.csv(celltype_fine,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_HumanPrimaryCellAtlasData_celltype_singleR_fine.csv"),row.names = F)


## DatabaseImmuneCellExpressionData
#refdata2 <-  DatabaseImmuneCellExpressionData()
#saveRDS(refdata2,"~/onedrive/Work/phD/phd_project/Gender_imm/results/DatabaseImmuneCellExpressionData.rds")
refdata2 <- readRDS("/public/home/ac8q5y82z9/data/3.Results/cell_identity/DatabaseImmuneCellexpressionData.rds")


cellpred_main <- SingleR(test = testdata1, ref = refdata2, labels = refdata2$label.main, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype=cellpred_main$labels, stringsAsFactors = F)
write.csv(celltype_main,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_DatabaseImmuneCellExpressionData_celltype_singleR.csv",row.names = F))

cellpred_fine <- SingleR(test = testdata1, ref = refdata2, labels = refdata2$label.fine, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype=cellpred_fine$labels, stringsAsFactors = F)
write.csv(celltype_fine,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_DatabaseImmuneCellExpressionData_celltype_singleR_fine.csv",row.names = F))


## MonacoImmuneData
#refdata3 <- MonacoImmuneData() 
#saveRDS(refdata3,"~/onedrive/Work/phD/phd_project/Gender_imm/results/MonacoImmuneData.rds")
refdata3 <- readRDS("/public/home/ac8q5y82z9/data/3.Results/cell_identity/MonacoImmuneData.rds")

clusters <- sce_cluster@meta.data$seurat_clusters
testdata1 <- GetAssayData(sce_cluster, slot="data")
cellpred_main <- SingleR(test = testdata1, ref = refdata3, labels = refdata3$label.main, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype=cellpred_main$labels, stringsAsFactors = F)
write.csv(celltype_main,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_MonacoImmuneData_celltype_singleR.csv",row.names = F))

cellpred_fine <- SingleR(test = testdata1, ref = refdata3, labels = refdata3$label.fine, 
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype=cellpred_fine$labels, stringsAsFactors = F)
write.csv(celltype_fine,paste0("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/res_0.5/res_0.5_MonacoImmuneData_celltype_singleR_fine.csv",row.names = F))


########### 12. Annotation #######
sce_cluster <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/res_0.5_sce_remove_doublet_tsne_clustered.rds")
sce_anno <- sce_cluster
anno_table <- data.frame("cluster" = c(0,1:16),"celltype" = c("CD8+ Tem", "CD4+ Naive T", 
                                                              "CD8+ Tcm", "B cells cluster 2",
                                                              "CD4+ Treg","Plasma cells cluster 3",
                                                              "B cells cluster 1","Plasma cells cluster 4",
                                                              "Basophils","Monocytes/Macrophages", 
                                                              "Endothelial cells", "CD8+ CXCL13+ Tex",
                                                              "Plasma cells cluster 1","Fibroblast",
                                                              "Plasma cells cluster 2","Plasma cells cluster 5",
                                                              "CD4+ Tfh"))


sce_anno@meta.data$barcode =rownames(sce_anno@meta.data)
sce_anno@meta.data=merge(sce_anno@meta.data,anno_table,by.x="seurat_clusters",by.y="cluster")
rownames(sce_anno@meta.data)= sce_anno@meta.data$barcode


#sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/sce_anno_res_0.5/res_0.5.rds")
#pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/annotated_dimplot.pdf")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/celltyep_remove_doubletonly_imm_spatialDimPlot.pdf",width = 8)
DimPlot(sce_anno, reduction = "tsne", group.by = "celltype", 
        pt.size=0.5, raster = F, cols = celltype_cols)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
DimPlot(sce_anno, reduction = "tsne", group.by = "orig.ident", 
        pt.size=0.5, raster = F)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
dev.off()
saveRDS(sce_anno,"~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")
sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")

########### 13. All marker distribution (FeaturePlot ) ########
#### Feature plot ###
marker_list <- c("COL1A1", "COL1A2", "CD19", "MS4A1", "CD79B", "CD79A","CST3",
                 "LYZ","HPGDS","MS4A2","CD3E","MZB1")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/markers_feature_plot.pdf",width = 18, height = 10)
cols = c("gray", "coral2")
FeaturePlot(sce_anno,reduction = "tsne", features = marker_list,cols = cols, pt.size = 1,ncol = 4 )+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/CD74_test_markers_feature_plot.pdf",width = 18, height = 10)
cols = c("gray", "coral2")
FeaturePlot(sce_anno,reduction = "tsne", features = , pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
dev.off()

pdf("/public/home/ac8q5y82z9/data/3.Results/SiT/merged/QC1/remove_doublet/cell_markers_Vlnplot_qc1.pdf",width = 15, height = 15)
StackedVlnPlot(sce_anno, marker_list, pt.size=0)
dev.off()

########### 14. ratio ###########
celltype_cols <- c( "#51c4c2", "#0d8a8c","#4583b3",
                      "#c63596" , "#be86ba", "#8b66b8",
                      "#f78e26",
                      "#f172ad", "#f7afb9",
                      "#4068b2","#512a93","#223271"
)
celltype_cols <- colorRampPalette(celltype_cols)(17)
scales::show_col(alpha(celltype_cols,1))
names(celltype_cols) <- levels(sce_anno$celltype)

Cellratio <- as.data.frame(prop.table(table(sce_anno$celltype,sce_anno$orig.ident), margin = 2))#计算各组样本不同细胞群比例

piechart_treatment <- function(Cellratio, treatment,col){
  cur_df <- Cellratio[Cellratio$Var2 == treatment,]
  cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
  pie = ggplot(cur_df, aes(x="", y=Freq, fill=Var1)) + geom_bar(stat="identity", width=1)+
    scale_fill_manual(values = col)
  
  
  # Convert to pie (polar coordinates) and add labels
  pie = pie + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))
  
  # Remove labels and add title
  pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = treatment)
  
  # Tidy up the theme
  pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank(),
                                      plot.title = element_text(hjust = 0.5, color = "#666666"))
  return(pie)
}

pie_before <- piechart_treatment(Cellratio, "before", celltype_cols)
pie_after <- piechart_treatment(Cellratio, "after", celltype_cols)
pie_after_ln <- piechart_treatment(Cellratio, "after_ln", celltype_cols)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/ratio_celltype_piechart.pdf", width = 8, height = 8)
ggarrange(pie_before, pie_after, pie_after_ln,common.legend = T,ncol = 3)
dev.off()

### sample vs celltype
cell_ratio_sample <- prop.table(table(sce_anno_tam_labeled$orig.ident, sce_anno_tam_labeled$celltype_combined), margin = 2)
cell_ratio_sample <- as.data.frame(cell_ratio_sample)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cell_ratio.pdf",width = 5, height = 5)
ggplot(cell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_jco()
  #theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) + 
dev.off()

### sample vs celltype in imm
sce_anno_tcell <- sce_anno@meta.data[sce_anno@meta.data$celltype %in% c("CD4+ Naive T", 
                                                                        "CD4+ Treg",
                                                                        "CD4+ Tfh","CD8+ Tem",  
                                                                        "CD8+ Tcm", 
                                                                        "CD8+ Effctor T"),]
cell_ratio_tcell <- prop.table(table(sce_anno_tcell$celltype,sce_anno_tcell$orig.ident), margin = 2)
cell_ratio_tcell <- as.data.frame(cell_ratio_tcell)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/cell_ratio_tcell.pdf",width = 5, height = 5)
ggplot(cell_ratio_tcell) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  #scale_fill_viridis(discrete = 1,option = "B")
  scale_fill_jco()
   theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) 
dev.off()  

ggbarstats(sce_anno@meta.data, celltype.y, orig.ident, palette = 'Set2',ggstatsplot.layer = FALSE)


 
##### T cells
sce_anno$t_cells <- ifelse(sce_anno@meta.data$celltype %in% c("CD4+ Naive T", 
                                                          "CD4+ Treg",
                                                          "CD4+ Tfh","CD8+ Tem",  
                                                          "CD8+ Tcm", 
                                                          "CD8+ Effctor T"),"T cells", "Non T cells")
tcell_ratio_sample <- prop.table(table(sce_anno$t_cells,sce_anno$orig.ident), margin = 2)
tcell_ratio_sample <- as.data.frame(tcell_ratio_sample)
tcell_ratio_sample$Var2 <- factor(tcell_ratio_sample$Var2, levels = c('before',"after", "after_ln"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_ratio.pdf",width = 5, height = 5)
ggplot(tcell_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_jco()

ggplot(tcell_ratio_sample[tcell_ratio_sample$Var1 == "T cells", ], aes(x =Var2, y= Freq, group = Var1)) + geom_line()
#theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) + 
dev.off()


##### Only CD8 CXCL13 
sce_anno$cd8_effector <- ifelse(sce_anno@meta.data$celltype == "CD8+ Effctor T" ,"CD8+ Effector T cells", "Other cells")
cd8_effector_ratio_sample <- prop.table(table(sce_anno$cd8_effector,sce_anno$orig.ident), margin = 2)
cd8_effector_ratio_sample <- as.data.frame(cd8_effector_ratio_sample)
cd8_effector_ratio_sample$Var2 <- factor(cd8_effector_ratio_sample$Var2, levels = c('before',"after", "after_ln"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cd8_effector_tcell_ratio.pdf",width = 5, height = 5)
ggplot(cd8_effector_ratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_jco()

ggplot(cd8_effector_ratio_sample[cd8_effector_ratio_sample$Var1 == "CD8+ Effector T cells", ],
       aes(x =Var2, y= Freq, group = Var1)) + geom_line(linetype="dotted", color="black", size= 1) +
  geom_point(size =5 ,color = "#00468B99") + ylab("Percentage of CD8+ effector T cells")
 
#theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) + 
dev.off()

########### 15. Hallmark functional analysis ####
## tcell 
tcells_only_sce_anno <- subset(sce_anno, subset = celltype %in% c("CD4+ Naive T", 
                                                                  "CD4+ Treg",
                                                                  "CD4+ Tfh",
                                                                  "CD8+ Tem",  
                                                                  "CD8+ Tcm", 
                                                                  "CD8+ Effctor T"))
before_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "before")
after_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "after")
after_ln_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "after_ln")
saveRDS(tcells_only_sce_anno, "~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcells_only_sce_anno.rds")

##### 15.1 all celltypes  #####
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_genesets <- subset(hallmark_genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
hallmark_genesets <- split(hallmark_genesets$gene_symbol, hallmark_genesets$gs_name)
## total
expr <- AverageExpression(sce_anno, assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)
gsva.res <- gsva(expr, hallmark_genesets, method = "gsva")
rownames(gsva.res) <- gsub("HALLMARK_","",rownames(gsva.res))
colnames(gsva.res) <- c("CD8+ Tem", "CD4+ Naive T", 
                        "CD8+ Tcm", "B cells unknown",
                        "CD4+ Treg","after unique 2",
                        "B cells","after_ln unique 1",
                        "Basophils","Monocytes", 
                        "Endothelial cells", "CD8+ Effctor T",
                        "Plasma cells cluster 1","Fibroblast",
                        "Plasma cells cluster 2","after_ln unique 2",
                        "CD4+ Tfh")

annotation_row <-  as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/hallmark.csv",col.names = F))
colnames(annotation_row) <- "Pathway"

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/tcell_heatmap_hallmark.pdf",width = 10, height = 10)
pheatmap(gsva.res[rownames(annotation_row),c("CD4+ Naive T", 
                                             "CD4+ Treg",
                                             "CD4+ Tfh",
                                             "CD8+ Tem",  
                                             "CD8+ Tcm", 
                                             "CD8+ Effctor T")],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()

#### test if CD74 coexpression with CD8 and CXCL13
colnames(expr) <- c("CD8+ Tem", "CD4+ Naive T", 
                        "CD8+ Tcm", "B cells unknown",
                        "CD4+ Treg","after unique 2",
                        "B cells","after_ln unique 1",
                        "Basophils","Monocytes", 
                        "Endothelial cells", "CD8+ Effctor T",
                        "Plasma cells cluster 1","Fibroblast",
                        "Plasma cells cluster 2","after_ln unique 2",
                        "CD4+ Tfh")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/CD74_test_markers_feature_heatmap.pdf",width = 18, height = 10)
pheatmap(expr[c(c("CD3D","CD3E","CD3G","CD8A","CD4","CXCL13","CD74","CD19","CD79A","CD79B","MS4A1")),],show_colnames = T, show_rownames = T, 
         scale = "row",
         angle_col = "45", cluster_rows = F,cluster_cols = T,
         #color = viridis(50,option = "D"),
         #color = col_fun,
         #annotation_row = annotation_row,
         #annotation_colors  =  annotation_colors,
         annotation_names_row = F) 

dev.off()

### sep

# Before
expr_before <- AverageExpression(before_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_before <- as.matrix(expr_before$RNA)
gsva.res_before <- gsva(expr_before, hallmark_genesets, method = "gsva")
rownames(gsva.res_before) <- gsub("HALLMARK_","",rownames(gsva.res_before))
colnames(gsva.res_before) <- c("CD8+ Tem", "CD4+ Naive T", 
                               "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                               "CD4+ Tfh")
gsva.res_before <- gsva.res_before[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                      "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]

expr_after <- AverageExpression(after_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after <- as.matrix(expr_after$RNA)
gsva.res_after <- gsva(expr_after, hallmark_genesets, method = "gsva")
rownames(gsva.res_after) <- gsub("HALLMARK_","",rownames(gsva.res_after))
colnames(gsva.res_after) <- c("CD8+ Tem", "CD4+ Naive T", 
                              "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                              "CD4+ Tfh")
gsva.res_after <- gsva.res_after[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                    "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]

expr_after_ln <- AverageExpression(after_ln_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after_ln <- as.matrix(expr_after_ln$RNA)
gsva.res_after_ln <- gsva(expr_after_ln, hallmark_genesets, method = "gsva")
rownames(gsva.res_after_ln) <- gsub("HALLMARK_","",rownames(gsva.res_after_ln))
colnames(gsva.res_after_ln) <- c("CD8+ Tem", "CD4+ Naive T", 
                                 "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                                 "CD4+ Tfh")
gsva.res_after_ln <- gsva.res_after_ln[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                          "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]



##### 15.2 T cell related only #####
### hallmark
annotation_row <-  as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/hallmark.csv",col.names = F))
colnames(annotation_row) <- "Pathway"
# total
expr <- AverageExpression(before_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)
gsva.res <- gsva(expr, hallmark_genesets, method = "gsva")
rownames(gsva.res) <- gsub("HALLMARK_","",rownames(gsva.res))
colnames(gsva.res) <- c("CD8+ Tem", "CD4+ Naive T", 
                               "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                               "CD4+ Tfh")

# Before
expr_before <- AverageExpression(before_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_before <- as.matrix(expr_before$RNA)
gsva.res_before <- gsva(expr_before, hallmark_genesets, method = "gsva")
rownames(gsva.res_before) <- gsub("HALLMARK_","",rownames(gsva.res_before))
colnames(gsva.res_before) <- c("CD8+ Tem", "CD4+ Naive T", 
                        "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                        "CD4+ Tfh")
gsva.res_before <- gsva.res_before[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                      "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]

expr_after <- AverageExpression(after_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after <- as.matrix(expr_after$RNA)
gsva.res_after <- gsva(expr_after, hallmark_genesets, method = "gsva")
rownames(gsva.res_after) <- gsub("HALLMARK_","",rownames(gsva.res_after))
colnames(gsva.res_after) <- c("CD8+ Tem", "CD4+ Naive T", 
                               "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                               "CD4+ Tfh")
gsva.res_after <- gsva.res_after[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                      "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]

expr_after_ln <- AverageExpression(after_ln_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after_ln <- as.matrix(expr_after_ln$RNA)
gsva.res_after_ln <- gsva(expr_after_ln, hallmark_genesets, method = "gsva")
rownames(gsva.res_after_ln) <- gsub("HALLMARK_","",rownames(gsva.res_after_ln))
colnames(gsva.res_after_ln) <- c("CD8+ Tem", "CD4+ Naive T", 
                               "CD8+ Tcm",  "CD4+ Treg", "CD8+ Effctor T",
                               "CD4+ Tfh")
gsva.res_after_ln <- gsva.res_after_ln[,c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                      "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" )]

### total
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/tcell_heatmap_hallmark.pdf",width = 10, height = 10)
pheatmap(gsva.res[rownames(annotation_row),c("CD4+ Naive T", 
                                             "CD4+ Treg",
                                             "CD4+ Tfh",
                                             "CD8+ Tem",  
                                             "CD8+ Tcm", 
                                             "CD8+ Effctor T")],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()
### sep
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/sep_tcell_heatmap_hallmark.pdf",width = 10, height = 10)
pheatmap(gsva.res_before[rownames(annotation_row),c("CD4+ Naive T", 
                                             "CD4+ Treg",
                                             "CD4+ Tfh",
                                             "CD8+ Tem",  
                                             "CD8+ Tcm", 
                                             "CD8+ Effctor T")],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)
pheatmap(gsva.res_after[rownames(annotation_row),c("CD4+ Naive T", 
                                                    "CD4+ Treg",
                                                    "CD4+ Tfh",
                                                    "CD8+ Tem",  
                                                    "CD8+ Tcm", 
                                                    "CD8+ Effctor T")],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)
pheatmap(gsva.res_after_ln[rownames(annotation_row),c("CD4+ Naive T", 
                                                    "CD4+ Treg",
                                                    "CD4+ Tfh",
                                                    "CD8+ Tem",  
                                                    "CD8+ Tcm", 
                                                    "CD8+ Effctor T")],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = T,
         color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()

### sep2
annotation_col <- data.frame(celltype = rep(colnames(gsva.res_before),3) , sample = rep(c("before","after","after_ln"),each = 6))

colnames(gsva.res_before) <- paste0(colnames(gsva.res_before),"_before")
colnames(gsva.res_after) <- paste0(colnames(gsva.res_after),"_after")
colnames(gsva.res_after_ln) <- paste0(colnames(gsva.res_after_ln),"_after_ln")

gsva.res_combine <- cbind(gsva.res_before,gsva.res_after,gsva.res_after_ln)
rownames(annotation_col)  <- colnames(gsva.res_combine)
annotation_col$celltype <- factor(annotation_col$celltype, levels =c("CD4+ Tfh","CD8+ Effctor T","CD4+ Treg",
                                                                     "CD8+ Tcm","CD4+ Naive T","CD8+ Tem" ) )
annotation_col <- annotation_col[order(annotation_col$celltype),]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/combine_tcell_heatmap_hallmark.pdf",width = 10, height = 10)

pheatmap(gsva.res_combine[rownames(annotation_row),rownames(annotation_col)],show_colnames = T, scale = "row", 
         angle_col = "45", cluster_rows = F,cluster_cols = F,
         color = viridis(100,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_col = annotation_col,
         gaps_col = c(3,6,9,12,15,18),
         annotation_names_row = F)
dev.off()




########### 16. Tcell related functional analysis ################
tcell_marker <- read_csv("rawData/new_tcell_marker.csv")
expr <- AverageExpression(tcells_only_sce_anno , assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)
expr <- expr[tcell_marker$Gene,]
#expr <- expr[,c("0","1","2","4","11")]

expr <- scale_mat(expr,"row")
dim(expr)
colnames(expr) <- c("CD8+ Tem","CD4+ Naive T","CD8+ Tcm","CD4+ Treg","CD8+ Effctor T","CD4+ Tfh")
annotation_row <- data.frame(tcell_marker$Group)
rownames(annotation_row) <- tcell_marker$Gene
dim(annotation_row)

annotation_colors <- viridis(13,option = "D")
names(annotation_colors) <- unique(tcell_marker$Group)
## total heatmap
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/only_tcell_related_heatmap.pdf",width = 6, height = 8)
pheatmap(expr,show_colnames = T, show_rownames = T, 
         scale = "row",
         angle_col = "45", cluster_rows = F,cluster_cols = T,
         #color = viridis(50,option = "D"),
         #color = col_fun,
         annotation_row = annotation_row,
         #annotation_colors  =  annotation_colors,
         annotation_names_row = F) 

dev.off()
### separate heatmap
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/sepapearate_only_tcell_related_heatmap.pdf",width = 6, height = 8)
pheatmap(expr,show_colnames = T, show_rownames = T, 
         scale = "row",
         angle_col = "45", cluster_rows = F,cluster_cols = T,
         #color = viridis(50,option = "D"),
         #color = col_fun,
         annotation_row = annotation_row,
         #annotation_colors  =  annotation_colors,
         annotation_names_row = F) 

dev.off()
### total Dotplot

split_tcell_marker <- split(tcell_marker$Gene,tcell_marker$Group)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/only_tcell_feature_dotplot.pdf",width = 12, height = 5)

DotPlot(tcells_only_sce_anno,features = split_tcell_marker[c("T","CD4","CD8", 
                                                             "Cytotoxic", "Exhausted","Naive",  "Tfh", "Treg")],
        group.by = "celltype",dot.scale = 7,scale.min = 10,scale = F) +
              theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) + scale_color_viridis(option = "B")
dev.off()
### split by sample
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/separate_only_tcell_feature_dotplot.pdf",width = 12, height = 5)

DotPlot(before_tcells_only_sce_anno,features = split_tcell_marker[c("T","CD4","CD8", 
                                                             "Cytotoxic", "Exhausted","Naive",  "Tfh", "Treg")],
        group.by = "celltype",dot.scale = 7,scale.min = 10,scale = F) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) + scale_color_viridis(option = "B")
DotPlot(after_tcells_only_sce_anno,features = split_tcell_marker[c("T","CD4","CD8", 
                                                                    "Cytotoxic", "Exhausted","Naive",  "Tfh", "Treg")],
        group.by = "celltype",dot.scale = 7,scale.min = 10,scale = F) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) + scale_color_viridis(option = "B")
DotPlot(after_ln_tcells_only_sce_anno,features = split_tcell_marker[c("T","CD4","CD8", 
                                                                    "Cytotoxic", "Exhausted","Naive",  "Tfh", "Treg")],
        group.by = "celltype",dot.scale = 7,scale.min = 10,scale = F) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) + scale_color_viridis(option = "B")


dev.off()


#### test if CD74 only express in cd8_ cxcl13
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/Tcell_CD74_test_markers_feature_heatmap.pdf",width = 18, height = 10)
pheatmap(expr[c("CD3D","CD3E","CD3G","CD8A","CD4","CXCL13","CD74"),],show_colnames = T, show_rownames = T, 
         scale = "row",
         angle_col = "45", cluster_rows = F,cluster_cols = T,
         #color = viridis(50,option = "D"),
         #color = col_fun,
         #annotation_row = annotation_row,
         #annotation_colors  =  annotation_colors,
         annotation_names_row = F) 

dev.off()



### check cluster 4 vs cluster 15
cluster_diff_markers <- FindMarkers(tcells_only_sce_anno,ident.1 = "4", ident.2 = "16")
cluster_diff_markers <- cluster_diff_markers[cluster_diff_markers$p_val <= 0.05,]
cluster_diff_markers <- cluster_diff_markers[cluster_diff_markers$avg_log2FC > 0,]

cluster_diff_markers[order(cluster_diff_markers$avg_log2FC,decreasing = T),]
dim(cluster_diff_markers)
cluster_diff_markers %>% top_n(n = 10, wt = abs(avg_log2FC))
########### 17. Correlation analysis #######
### All celltypes
expr <- AverageExpression(sce_anno,group.by = "celltype", assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

top1000_genes <- names(tail(sort(apply(expr,1,sd)),1000))
#View(expr[top1000_genes,])
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/correlation_all_celltypes_heatmap.pdf",width = 6, height = 8)

pheatmap(cor(expr[top1000_genes,], method = "spearman"))
corrplot::corrplot(cor(expr[top1000_genes,], method = "spearman"), method = "pie")

dev.off()

### only Tcells (all samples together)
expr <- AverageExpression(tcells_only_sce_anno ,group.by = "celltype", assays = "RNA", slot = "data")
expr <- as.matrix(expr$RNA)

top1000_genes <- names(tail(sort(apply(expr,1,sd)),1000))
#View(expr[top1000_genes,])
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/correlation_Tcells_heatmap.pdf",width = 6, height = 8)

pheatmap(cor(expr[top1000_genes,], method = "spearman"))
corrplot::corrplot(cor(expr[top1000_genes,], method = "spearman"), method = "pie")

dev.off()

### only T cells (sep)
expr_before <- AverageExpression(before_tcells_only_sce_anno, group.by = "celltype", assays = "RNA", slot = "data")
expr_before <- as.matrix(expr_before$RNA)

expr_after <- AverageExpression(after_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after <- as.matrix(expr_after$RNA)

expr_after_ln <- AverageExpression(after_ln_tcells_only_sce_anno, assays = "RNA", slot = "data")
expr_after_ln <- as.matrix(expr_after_ln$RNA)

annotation_col <- data.frame(celltype = rep(colnames(expr_before),3) , sample = rep(c("before","after","after_ln"),each = 6))

# annotation                                        
colnames(expr_before) <- paste0(colnames(expr_before),"_before")
colnames(expr_after) <- paste0(colnames(expr_after),"_after")
colnames(expr_after_ln) <- paste0(colnames(expr_after_ln),"_after_ln")

expr_combine <- cbind(expr_before,expr_after,expr_after_ln)
cor_matrix <- cor(expr_combine[top1000_genes,], method = "spearman")

rownames(annotation_col)  <- colnames(expr_combine)


pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/sample_combine_tcell_heatmap_hallmark.pdf",width = 12, height = 10)
pheatmap(cor_matrix,annotation_row = annotation_col,annotation_col = annotation_col,angle_col = "45") 
dev.off()


