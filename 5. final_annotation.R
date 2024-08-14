
before@meta.data$celltype <- ifelse(before$seurat_clusters == "0",'Epithelial cells (KRT17+)',
                                    ifelse(before$seurat_clusters == "1",'Fibroblast',
                                           ifelse(before$seurat_clusters == "2",'Epithelial cells (CSTB+), T cells, B cells, Macrophage',
                                                  ifelse(before$seurat_clusters == "3",'Bcells T cells',
                                                         ifelse(before$seurat_clusters == "4",'Malignant Epithelial cells',
                                                                ifelse(before$seurat_clusters == "5",'Epithelial cells (KRT17+, CSTB+), Macrophage',"no"))))))


after@meta.data$celltype <- ifelse(after$seurat_clusters == "0",'Unknown',
                                   ifelse(after$seurat_clusters == "1",'Fibroblast',
                                          ifelse(after$seurat_clusters == "2",'Tcell Macrophage Fibroblast ',
                                                 ifelse(after$seurat_clusters == "3",'Epithelial cells // Unknown',
                                                        ifelse(after$seurat_clusters == "4",'Epithelial cells (CSTB+) and Macrophage',
                                                               ifelse(after$seurat_clusters == "5",'Malignant Epithelial cells , Macrophage, T cell',
                                                                      ifelse(after$seurat_clusters == "6",'B cells, Tcells, Fibroblast,Macrophage',
                                                                             ifelse(after$seurat_clusters == "7",'Epithelial cells (KRT17+,CSTB+) , Macrophage',"no"))))))))

after_ln@meta.data$celltype <- ifelse(after_ln$seurat_clusters == "0","Tcells and Macrophage cluster 1",
                                      ifelse(after_ln$seurat_clusters == "1","Fibroblast, Tcells and Macrophage",
                                             ifelse(after_ln$seurat_clusters == "2","Tcells and Macrophage cluster 2",
                                                    ifelse(after_ln$seurat_clusters == "3","Tcells and Macrophage cluster 3",
                                                           ifelse(after_ln$seurat_clusters == "4","Epithelial cells and Tcells, Macrophage",
                                                                  ifelse(after_ln$seurat_clusters == "5","Malignant Epithelial cells, Tcells, Macrophage","oh no"))))))

###  Dimplot (Fig 3A) ####
pdf("~/onedrive/Work/phD/phd_project/SiT/results/celltype_spatial_dimplot.pdf")
SpatialDimPlot(after_ln,group.by = "celltype", pt.size.factor=1.2, alpha = c(1, 1))
SpatialDimPlot(after,group.by = "celltype", pt.size.factor=1.8, alpha = c(1, 1)) 
SpatialDimPlot(before,group.by = "celltype", pt.size.factor=2.9, alpha = c(1, 1)) 
dev.off()

