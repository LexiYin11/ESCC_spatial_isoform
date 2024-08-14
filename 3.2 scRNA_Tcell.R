#### 0. Library import #####
library(monocle3)
#### 1. Extract  Tcell only ######
sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")

tcells_only_sce_anno_ori <- subset(sce_anno, subset = celltype %in% c("CD4+ Naive T", 
                                                                  "CD4+ Treg",
                                                                  "CD4+ Tfh",
                                                                  "CD8+ Tem",  
                                                                  "CD8+ Tcm", 
                                                                  "CD8+ Effctor T"))
before_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "before")
after_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "after")
after_ln_tcells_only_sce_anno <- subset(tcells_only_sce_anno, subset = orig.ident == "after_ln")


tcells_only_sce_anno <- SCTransform(tcells_only_sce_anno_ori) %>% RunPCA() 
tcells_only_sce_anno <- RunHarmony(tcells_only_sce_anno, group.by.vars = "orig.ident",
                                   assay.use = "SCT", max.iter.harmony = 10)  
## 降维聚类
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
tcells_only_sce_anno <- RunUMAP(tcells_only_sce_anno, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.8)   
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_dimplot.pdf",width = 8)

DimPlot(tcells_only_sce_anno, reduction = "umap", group.by = "celltype", label = F) 
DimPlot(tcells_only_sce_anno, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(tcells_only_sce_anno, reduction = "umap", group.by = "seurat_clusters", label = F) 

DimPlot(tcells_only_sce_anno_ori, reduction = "umap", group.by = "celltype", label = F) 
DimPlot(tcells_only_sce_anno_ori, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(tcells_only_sce_anno_ori, reduction = "umap", group.by = "seurat_clusters", label = F) 


dev.off()
saveRDS(tcells_only_sce_anno_ori,"~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcells_only_sce_anno_ori.rds")
#### 2. Monocle3 #######
### 2.1 total ####
tcells_only_sce_anno <- tcells_only_sce_anno_ori
data <- GetAssayData(tcells_only_sce_anno, assay = 'RNA', slot = 'counts')
cell_metadata <- tcells_only_sce_anno@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(tcells_only_sce_anno, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
## 
cds <- learn_graph(cds)
p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE)
p4 = plot_cells(cds, color_cells_by = "celltype",label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE)

p = wrap_plots(p1, p2,p3)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_monocle3.pdf",width = 8)
p1
p2
p3
p4
dev.off()

### choose root

cds <- order_cells(cds)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_monocle3_Trajectory_Pseudotime.pdf",width = 8)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()

### save rds 
saveRDS(cds, "~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cds.rds" )
### differential analysis
track_genes <- graph_test(cds, neighbor_graph = "principal_graph",cores = 6)
track_genes <- track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

write.csv(track_genes,"~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_monocle3__Trajectory_Pseudotime.csv")
track_genes <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_monocle3__Trajectory_Pseudotime.csv")
### top10
track_genes_sig <- track_genes %>% top_n(n=30, morans_I) %>% pull(gene_short_name) %>% as.character()
genes <- c("CD8A", "CXCL13", "CD74", "CCL4","CCL5", "GZMK")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcell_monocle3_Trajectory_Pseudotime_genes_jitterplot.pdf",width = 8,height = 4)
plot_genes_in_pseudotime(cds[genes,], color_cells_by = "celltype",
                         min_expr=0.5, ncol=3)
dev.off()

### add pseudotime label
pseudotime <- pseudotime(cds, reduction_method = "UMAP")
pseudotime <- pseudotime[rownames(tcells_only_sce_anno@meta.data)]
tcells_only_sce_anno$psedotime <- pseudotime

### 2.2 Sep ####
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/sample_tcell_dimplot.pdf",width = 8)

DimPlot(before_tcells_only_sce_anno, reduction = "umap", group.by = "celltype", label = F) 
#DimPlot(before_tcells_only_sce_anno, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(before_tcells_only_sce_anno, reduction = "umap", group.by = "seurat_clusters", label = F) 

DimPlot(after_tcells_only_sce_anno, reduction = "umap", group.by = "celltype", label = F) 
#DimPlot(after_tcells_only_sce_anno, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(after_tcells_only_sce_anno, reduction = "umap", group.by = "seurat_clusters", label = F) 

DimPlot(after_ln_tcells_only_sce_anno, reduction = "umap", group.by = "celltype", label = F) 
#DimPlot(after_ln_tcells_only_sce_anno, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(after_ln_tcells_only_sce_anno, reduction = "umap", group.by = "seurat_clusters", label = F) 

dev.off()
## Before
##创建CDS对象并预处理数据
tcells_only_sce_anno <- before_tcells_only_sce_anno
data <- GetAssayData(tcells_only_sce_anno, assay = 'RNA', slot = 'counts')
cell_metadata <- tcells_only_sce_anno@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(tcells_only_sce_anno, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
## 识别轨迹
cds <- learn_graph(cds)
p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE)

#p3 <- p3 + geom_vline(xintercept = seq(-5,-4,0.25)) + geom_hline(yintercept = seq(0,2,0.25))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/before_tcell_monocle3.pdf",width = 8)
p1
p2
p3
dev.off()

embed <- data.frame(Embeddings(tcells_only_sce_anno, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -4.75 & UMAP_1 < -4.5 & UMAP_2 > 0.25 & UMAP_2 < 0.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)

## after
##创建CDS对象并预处理数据
tcells_only_sce_anno <- after_tcells_only_sce_anno
data <- GetAssayData(tcells_only_sce_anno, assay = 'RNA', slot = 'counts')
cell_metadata <- tcells_only_sce_anno@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(tcells_only_sce_anno, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
## 识别轨迹
cds <- learn_graph(cds)
p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE)

#p3 <- p3 + geom_vline(xintercept = seq(-5,-4,0.25)) + geom_hline(yintercept = seq(0,2,0.25))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/after_tcell_monocle3.pdf",width = 8)
p1
p2
p3
dev.off()

embed <- data.frame(Embeddings(tcells_only_sce_anno, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -4.75 & UMAP_1 < -4.5 & UMAP_2 > 0.25 & UMAP_2 < 0.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)

## after_ln
##创建CDS对象并预处理数据
tcells_only_sce_anno <- after_ln_tcells_only_sce_anno
data <- GetAssayData(tcells_only_sce_anno, assay = 'RNA', slot = 'counts')
cell_metadata <- tcells_only_sce_anno@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(tcells_only_sce_anno, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
## 识别轨迹
cds <- learn_graph(cds)
p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                label_branch_points = FALSE)

#p3 <- p3 + geom_vline(xintercept = seq(-5,-4,0.25)) + geom_hline(yintercept = seq(0,2,0.25))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/after_ln_tcell_monocle3.pdf",width = 8)
p1
p2
p3
dev.off()

embed <- data.frame(Embeddings(tcells_only_sce_anno, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -4.75 & UMAP_1 < -4.5 & UMAP_2 > 0.25 & UMAP_2 < 0.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)


#### 3.Extract CD8 effector only #######
cd8_tcells_only_sce_anno_ori <- subset(sce_anno, subset = celltype %in% c("CD8+ Effctor T"))
cd8_tcells_only_sce_anno <- SCTransform(cd8_tcells_only_sce_anno_ori) %>% RunPCA() 
cd8_tcells_only_sce_anno <- RunHarmony(cd8_tcells_only_sce_anno, group.by.vars = "orig.ident",
                                   assay.use = "SCT", max.iter.harmony = 10)  
## 降维聚类
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
cd8_tcells_only_sce_anno <- RunUMAP(cd8_tcells_only_sce_anno, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.8)   
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cd8_effector_tcell_dimplot.pdf",width = 8)

DimPlot(cd8_tcells_only_sce_anno, reduction = "umap", group.by = "celltype", label = F) 
DimPlot(cd8_tcells_only_sce_anno, reduction = "umap", group.by = "orig.ident", label = F) 
DimPlot(cd8_tcells_only_sce_anno, reduction = "umap", group.by = "seurat_clusters", label = F) 

dev.off()

#### 4.DEG analysis #####

Idents(cd8_tcells_only_sce_anno) = cd8_tcells_only_sce_anno$orig.ident
table(Idents(cd8_tcells_only_sce_anno))
### before vs after vs after_ln 
all_markers <- FindAllMarkers(cd8_tcells_only_sce_anno)
                              #[,cd8_tcells_only_sce_anno$orig.ident==x],ident.1 = 'c1',
              #ident.2 = 'c2')
write.csv(all_markers,"~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cd8_effector_markers.csv")
top50 <- all_markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
write.csv(top50,"~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cd8_effector_markers_top50.csv")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/cd8_effector_markers_heatmap.pdf",width = 8)

DoHeatmap(cd8_tcells_only_sce_anno ,top10$gene)
dev.off()
