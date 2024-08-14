####### 0. Library import #######
library(SPATA2)
library(Seurat)
library(ggsci)
library(viridis)
library(msigdbr)
library(GSVA)
library(matrixStats)
library(ggstatsplot)
library(ggpubr)
library(dplyr)

####### 1. Data import #######
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_ratio_labeled.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_ratio_labeled.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")

after_ln_iso_cluster <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_isoform_clustering.rds")
after_iso_cluster <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_isoform_clustering.rds")
before_iso_cluster <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_isoform_clustering.rds")

####### 2. Quantify the contact frequency (Fig 3B) #######
# Each cell has its own contact frequency and calculate the mean value 
##### 2.1 Functions
enrichment_function_list_score <- function(seu_obj, cell_list,  coordinates){
  enrichment_label_list <- c()
  metadata <- seu_obj@meta.data
  for (cell in cell_list){
    #label <- metadata[cell,type]
    enrichment_label_list <- c(enrichment_label_list, 
                               enrichment_function_score(seu_obj, cell,  coordinates))
  }
  return(enrichment_label_list)
}
enrichment_function_score <- function(seu_obj, cell, coordinates){
  metadata <- seu_obj@meta.data
  cell_cluster <- metadata[metadata$cell == cell, "seurat_clusters"]
  cell_neighbors<- get_neighbor_function(coordinates,cell=cell)
  
  cell_neighbors_meta <- metadata[metadata$cell %in% cell_neighbors,]
  return(sum(cell_neighbors_meta$seurat_clusters != cell_cluster))
    
}
get_neighbor_function <- function(coordinates, cell){
  neighbors <- c()
  x <- coordinates[coordinates$barcode == cell,"x"]
  y <- coordinates[coordinates$barcode == cell,"y"]
  neighbors_coordinate <- list(c(x - 2 , y), c(x + 2 , y),
                               c(x - 1 , y + 1), c(x +1, y+ 1),
                               c(x -1 , y-1),  c(x+1, y- 1))
  #twice_neighbors <- c()
  # for (xy in neighbors_coordinate){
  #   x1 = xy[[1]]
  #   y1 = xy[[2]]
  #   neighbors_ <- list(c(x1 - 2 , y1), c(x1 + 2 , y1),
  #                     c(x1 - 1 , y1 + 1), c(x1 +1, y1+ 1),
  #                     c(x1 -1 , y1-1),  c(x1+1, y1- 1))
  #   twice_neighbors <- c(twice_neighbors, neighbors_)
  # }
  for (xy in neighbors_coordinate){
    coordinates_df <- coordinates %>% filter(x == xy[1], y == xy[2])
    if (nrow(coordinates_df) !=0){
      neighbors <- c(neighbors, coordinates_df$barcode)
    }
  }
  return(unique(neighbors))
}
##### 2.2 after_ln
after_ln_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after_LN/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  dplyr::select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
## score frequency
after_ln$contact_frequency <- enrichment_function_list_score(after_ln,after_ln$cell,  after_ln_coordinates)
summary(after_ln$contact_frequency)
after_ln_contact <- after_ln@meta.data[,c("cell", "contact_frequency")]
after_ln_contact$sample <- "after_ln"
##### 2.3 after
after_coordinates <-  LoadSpatialCoordinats("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_After/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  dplyr::select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))

## score frequency
after$contact_frequency <- enrichment_function_list_score(after,after$cell,  after_coordinates)
summary(after$contact_frequency)
after_contact <- after@meta.data[,c("cell", "contact_frequency")]
after_contact$sample <- "after"
##### 2.4 before
before$cell <- rownames(before@meta.data)
before_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/1.Count/P1_before/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  dplyr::select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))

## score frequency
before$contact_frequency <- enrichment_function_list_score(before,before$cell,  before_coordinates)
summary(before$contact_frequency)
before_contact <- before@meta.data[,c("cell", "contact_frequency")]
before_contact$sample <- "before"
##### 2.5 Visualization 
freq_contact <- rbind(before_contact, after_contact)
freq_contact <- rbind(freq_contact, after_ln_contact)
freq_contact$sample <- factor(freq_contact$sample, levels = c("before", "after", "after_ln"))
p1 <- ggplot(freq_contact , aes(x= `sample`, y=`contact_frequency`, fill=factor(`sample`)))+ 
  geom_boxplot(outlier.shape=NA) + 
  ylab("Contact frequency") +
  xlab("Sample") + theme_classic() + theme(legend.position = "none") +  
  stat_compare_means(method = "t.test",comparisons = list(c("before", "after"), c("before", "after_ln"),
                                                          c("after", "after_ln")))
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation_contact_freq.pdf")
p1
dev.off()
saveRDS(before, "~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_ratio_labeled.rds")
saveRDS(after, "~/onedrive/Work/phD/phd_project/SiT/results/before_minion/after_ratio_labeled.rds")
saveRDS(after_ln, "~/onedrive/Work/phD/phd_project/SiT/results/before_minion/after_ln_ratio_labeled.rds")


####### 3. Calculation of similarity of isoform and gene clustering (Fig S2) #######

##### 3.1 Functions
enrichment_function_list_similarity <- function(seu_obj, cell_list,  coordinates){
  enrichment_label_list <- c()
  metadata <- seu_obj@meta.data
  for (cell in cell_list){
    #label <- metadata[cell,type]
    enrichment_label_list <- c(enrichment_label_list, 
                               enrichment_function_similarity(seu_obj, cell,  coordinates))
  }
  return(enrichment_label_list)
}
enrichment_function_similarity <- function(seu_obj, cell, coordinates){
  metadata <- seu_obj@meta.data
  cell_metadata <- metadata[metadata$cell == cell,]
  cell_iso_cluster <- cell_metadata$ISO_clusters
  cell_neighbors <- get_neighbor_function(coordinates,cell=cell)
  total_metadata <- rbind(cell_metadata, metadata[metadata$cell %in% cell_neighbors,])
  return(ifelse(cell_iso_cluster %in% total_metadata$SCT_snn_res.0.5, 1, 0 ))
  
}

##### 3.2 after_ln
## replace to same number
after_ln_iso_cluster$ISO_clusters <- after_ln_iso_cluster$ISO_snn_res.0.7
after_ln_iso_cluster$ISO_clusters <- ifelse(after_ln_iso_cluster$ISO_snn_res.0.7 == 1, 3,
                                            ifelse(after_ln_iso_cluster$ISO_snn_res.0.7 == 2, 1,
                                                   ifelse(after_ln_iso_cluster$ISO_snn_res.0.7 == 3, 2,
                                                          ifelse(after_ln_iso_cluster$ISO_snn_res.0.7 == 4, 5,
                                                                 ifelse(after_ln_iso_cluster$ISO_snn_res.0.7 == 5, 4, 0)))))

  
head(after_ln_iso_cluster@meta.data)
SpatialDimPlot(after_ln_iso_cluster, group.by = "ISO_clusters")
SpatialDimPlot(after_ln_iso_cluster, group.by = "SCT_snn_res.0.5")
## coordinates
after_ln_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after_LN/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  dplyr::select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
## similarity
after_ln_iso_cluster$cell <- rownames(after_ln_iso_cluster@meta.data)
after_ln_iso_cluster$next_third_similarity <- enrichment_function_list_similarity(after_ln_iso_cluster,after_ln_iso_cluster$cell,  after_ln_coordinates)
sum(after_ln_iso_cluster$next_third_similarity) / nrow(after_ln_iso_cluster@meta.data)
## 0.7632664

##### 3.3 after
## replace to same number
after_iso_cluster$ISO_clusters <- after_iso_cluster$ISO_snn_res.0.9
after_iso_cluster$ISO_clusters <- ifelse(after_iso_cluster$ISO_snn_res.0.9 == 1, 7,
                                            ifelse(after_iso_cluster$ISO_snn_res.0.9 == 2, 3,
                                                   ifelse(after_iso_cluster$ISO_snn_res.0.9 == 3, 5,
                                                          ifelse(after_iso_cluster$ISO_snn_res.0.9 == 4, 0,
                                                                 ifelse(after_iso_cluster$ISO_snn_res.0.9 == 5, 1, 
                                                                        ifelse(after_iso_cluster$ISO_snn_res.0.9 == 6, 6,
                                                                               ifelse(after_iso_cluster$ISO_snn_res.0.9 == 7, 2,
                                                                                      ifelse(after_iso_cluster$ISO_snn_res.0.9 == 8, 8,4))))))))


head(after_iso_cluster@meta.data)
SpatialDimPlot(after_iso_cluster, group.by = "ISO_clusters")
SpatialDimPlot(after_iso_cluster, group.by = "SCT_snn_res.0.5")

## coordinates
after_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  dplyr::select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
## similarity
after_iso_cluster$cell <- rownames(after_iso_cluster@meta.data)
after_iso_cluster$next_third_similarity <- enrichment_function_list_similarity(after_iso_cluster,after_iso_cluster$cell,  after_coordinates)
sum(after_iso_cluster$next_third_similarity) / nrow(after_iso_cluster@meta.data)
## 0.7299035

##### 3.4 Before
## replace to same number
before_iso_cluster$ISO_clusters <- before_iso_cluster$ISO_snn_res.0.8
before_iso_cluster$ISO_clusters <- ifelse(before_iso_cluster$ISO_snn_res.0.8 == 1, 5,
                                         ifelse(before_iso_cluster$ISO_snn_res.0.8 == 2, 2,
                                                ifelse(before_iso_cluster$ISO_snn_res.0.8 == 3, 0,
                                                       ifelse(before_iso_cluster$ISO_snn_res.0.8 == 4, 1,
                                                              ifelse(before_iso_cluster$ISO_snn_res.0.8 == 5, 4, 3)))))


head(before_iso_cluster@meta.data)
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation_similarity_spatial_dimplot.pdf")
SpatialDimPlot(before_iso_cluster, group.by = "ISO_clusters",pt.size.factor = 3)
SpatialDimPlot(before_iso_cluster, group.by = "SCT_snn_res.0.5",pt.size.factor = 3)
dev.off()
## similarity
before_iso_cluster$cell <- rownames(before_iso_cluster@meta.data)
before_iso_cluster$next_third_similarity <- enrichment_function_list_similarity(before_iso_cluster,before_iso_cluster$cell,  before_coordinates)
sum(before_iso_cluster$next_third_similarity) / nrow(before_iso_cluster@meta.data)

## 0.5873984


####### 4. Check gene expression exclusive (reviewer response) ######
##### 4.1 import single cell data 
tcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcells_only_sce_anno_ori.rds")
TAM_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/macrophage/old/res_0.5_TAM_only_sce_anno_anno.rds")

##### 4.2  T cells ####
cd8_cxcl13_tex <- subset(tcells_only_sce_anno, subset = celltype == "CD8+ Effctor T")
expr_tcell <- cd8_cxcl13_tex@assays$RNA@counts[c("CD74","C1QC"),]
expr_tcell <- as.matrix(expr_tcell)
expr_tcell <- as.data.frame(t(expr_tcell))
summary(expr_tcell$CD74)
summary(expr_tcell$C1QC)

label_tcell <- expr_tcell
label_tcell$CD74 <- ifelse(label_tcell$CD74 == 0, "Unexpressed", "Expressed")
label_tcell$C1QC <- ifelse(label_tcell$C1QC == 0, "Unexpressed", "Expressed")

### 4.2.1 CD74
ratio_cd74 <- prop.table(table(label_tcell$CD74) )#计算各组样本不同细胞群比例
cur_df <- as.data.frame(t(ratio_cd74))
cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
pie_cd74 = ggplot(cur_df, aes(x="", y=Freq, fill=Var2)) + geom_bar(stat="identity", width=1)+
  scale_fill_manual(values = c( "#370857","grey"))


# Convert to pie_cd74 (polar coordinates) and add labels
pie_cd74 = pie_cd74 + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))

# Remove labels and add title
pie_cd74 = pie_cd74 + labs(x = NULL, y = NULL, fill = NULL)

# Tidy up the theme
pie_cd74 = pie_cd74 + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "white"))
### 4.2.2 C1QC
ratio_c1qc <- prop.table(table(label_tcell$C1QC) )#计算各组样本不同细胞群比例
cur_df <- as.data.frame(t(ratio_c1qc))
cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
pie_c1qc = ggplot(cur_df, aes(x="", y=Freq, fill=Var2)) + geom_bar(stat="identity", width=1)+
  scale_fill_manual(values = c( "grey"))


# Convert to pie_c1qc (polar coordinates) and add labels
pie_c1qc = pie_c1qc + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))

# Remove labels and add title
pie_c1qc = pie_c1qc + labs(x = NULL, y = NULL, fill = NULL)

# Tidy up the theme
pie_c1qc = pie_c1qc + theme_classic() + theme(axis.line = element_blank(),
                                              axis.text = element_blank(),
                                              axis.ticks = element_blank(),
                                              plot.title = element_text(hjust = 0.5, color = "white"))





##### 4.3  TAM cells ####
c1qc_tam <- subset(TAM_only_sce_anno, subset = tam_label == "C1QC+ TAM")
c1qc_tam <- c1qc_tam@assays$RNA@counts[c("CXCL13","CD8A"),]
c1qc_tam <- as.matrix(c1qc_tam)
c1qc_tam <- as.data.frame(t(c1qc_tam))
summary(c1qc_tam$CXCL13)
summary(c1qc_tam$CD8A)

label_tcell <- c1qc_tam
label_tcell$CD8A <- ifelse(label_tcell$CD8A == 0, "Unexpressed", "Expressed")
label_tcell$CXCL13 <- ifelse(label_tcell$CXCL13 == 0, "Unexpressed", "Expressed")

### 4.3.1 CD8A
ratio_CD8A <- prop.table(table(label_tcell$CD8A) )#计算各组样本不同细胞群比例
cur_df <- as.data.frame(t(ratio_CD8A))
cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
pie_CD8A = ggplot(cur_df, aes(x="", y=Freq, fill=Var2)) + geom_bar(stat="identity", width=1)+
  scale_fill_manual(values = c( "#BB3754FF","grey"))


# Convert to pie_CD8A (polar coordinates) and add labels
pie_CD8A = pie_CD8A + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))

# Remove labels and add title
pie_CD8A = pie_CD8A + labs(x = NULL, y = NULL, fill = NULL)

# Tidy up the theme
pie_CD8A = pie_CD8A + theme_classic() + theme(axis.line = element_blank(),
                                              axis.text = element_blank(),
                                              axis.ticks = element_blank(),
                                              plot.title = element_text(hjust = 0.5, color = "white"))
### 4.3.2 CXCL13
ratio_cxcl13 <- prop.table(table(label_tcell$CXCL13) )#计算各组样本不同细胞群比例
cur_df <- as.data.frame(t(ratio_cxcl13))
cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
pie_cxcl13 = ggplot(cur_df, aes(x="", y=Freq, fill=Var2)) + geom_bar(stat="identity", width=1)+
  scale_fill_manual(values = c( "#BB3754FF","grey"))


# Convert to pie_cxcl13 (polar coordinates) and add labels
pie_cxcl13 = pie_cxcl13 + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))

# Remove labels and add title
pie_cxcl13 = pie_cxcl13 + labs(x = NULL, y = NULL, fill = NULL)

# Tidy up the theme
pie_cxcl13 = pie_cxcl13 + theme_classic() + theme(axis.line = element_blank(),
                                              axis.text = element_blank(),
                                              axis.ticks = element_blank(),
                                              plot.title = element_text(hjust = 0.5, color = "white"))

pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation/gene_exclusive_piechart.pdf")
ggarrange(pie_c1qc, pie_cd74,common.legend = T, ncol = 2)
ggarrange(pie_CD8A, pie_cxcl13, common.legend = T, ncol = 2)
dev.off()

##### 4.4  Visualization ####
CD74_tam <- c1qc_tam@assays$SCT@data["CD74",]
CD74_tam <- as.data.frame(CD74_tam)
cd8_cxcl13_tex_scaled <- ScaleData(cd8_cxcl13_tex,assay = "RNA")
CD74_tex <- cd8_cxcl13_tex_scaled@assays$RNA@data[c("CD74"),]
CD74_tex <- as.data.frame(CD74_tex)

CD74_tam$celltype <- "C1QC+ TAM"
colnames(CD74_tam)[1] <- "CD74"
CD74_tex$celltype <- "CD8+ CXCL13+ Tex"
colnames(CD74_tex)[1] <- "CD74"
cd74_total <- rbind(CD74_tam, CD74_tex)
## boxplot
p1 <- ggplot(cd74_total , aes(x= `celltype`, y=`CD74`, fill=factor(`celltype`)))+ 
  geom_boxplot(outlier.shape=NA) + 
  ylab("CD74 expression") +
  xlab("celltype") + theme_classic() + theme(legend.position = "none") +  
  stat_compare_means(method = "t.test")
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation/CD74_expression.pdf")
p1
dev.off()

####### 4.  RPS9 isoform 3D stack visualization (Fig 1E)######
isoforms_RPS9 <- c("RPS9..ENST00000391751","RPS9..ENST00000391753","RPS9..ENST00000302907")

## Before
DefaultAssay(before) <- "ISO"
p1 <- SpatialFeaturePlot(before, features = isoforms_RPS9[1], pt.size.factor = 3) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p2 <- SpatialFeaturePlot(before, features = isoforms_RPS9[2], pt.size.factor = 3) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p3 <- SpatialFeaturePlot(before, features = isoforms_RPS9[3], pt.size.factor = 3) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation/before_RPS9_expression.pdf")
ggarrange(p1,p2,p3, common.legend = T, ncol = 3)
dev.off()
## After
DefaultAssay(after) <- "ISO"
p1 <- SpatialFeaturePlot(after, features = isoforms_RPS9[1], pt.size.factor = 2) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p2 <- SpatialFeaturePlot(after, features = isoforms_RPS9[2], pt.size.factor = 2) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p3 <- SpatialFeaturePlot(after, features = isoforms_RPS9[3], pt.size.factor = 2) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation/after_RPS9_expression.pdf")
ggarrange(p1,p2,p3, common.legend = T, ncol = 3)
dev.off()
## After_ln
DefaultAssay(after_ln) <- "ISO"
p1 <- SpatialFeaturePlot(after_ln, features = isoforms_RPS9[1], pt.size.factor = 1.4) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p2 <- SpatialFeaturePlot(after_ln, features = isoforms_RPS9[2], pt.size.factor = 1.4) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
p3 <- SpatialFeaturePlot(after_ln, features = isoforms_RPS9[3], pt.size.factor = 1.4) + 
  scale_fill_viridis(option = "A",discrete = F,limits = c(0,4), begin = 0.2) 
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/2nd_validation/after_ln_RPS9_expression.pdf")
ggarrange(p1,p2,p3, common.legend = T, ncol = 3)
dev.off()