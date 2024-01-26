####### 0. Library prep ########
library(testSctpa)
library(msigdbr)
library(GSVA)
library(pheatmap)
####### 1. Data import #######
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_genesets <- subset(hallmark_genesets, select = c("gs_name", "gene_symbol")) %>% as.data.frame()
hallmark_genesets <- split(hallmark_genesets$gene_symbol, hallmark_genesets$gs_name)


####### 2. GSVA -- Hallmark #######
### before
expr_before <- AverageExpression(before, assays = "SCT", slot = "data")
expr_before <- as.matrix(expr_before$SCT)
gsva.res_before <- gsva(expr_before, hallmark_genesets, method = "gsva")
#write.csv(gsva.res, "~/onedrive/Work/phD/phd_project/SiT/results/after_ln/after_ln_hallmark_gsva.csv")
rownames(gsva.res_before) <- gsub("HALLMARK_","",rownames(gsva.res))
#rownames(gsva.res) <- gsub("_","", rownames(gsva.res))

colnames(gsva.res_before) <- c('Epithelial cells (KRT17+)','Fibroblast', 'Epithelial cells (CSTB+), T cells, B cells, Macrophage',
                        'Bcells T cells','Malignant Epithelial cells', 'Epithelial cells (KRT17+, CSTB+), Macrophage')
gsva.res_before <- gsva.res_before[,c('Fibroblast','Epithelial cells (KRT17+)', 
                        'Epithelial cells (KRT17+, CSTB+), Macrophage',
                        'Epithelial cells (CSTB+), T cells, B cells, Macrophage',
                        'Bcells T cells',
                        'Malignant Epithelial cells'
                        )]
### after
expr_after <- AverageExpression(after, assays = "SCT", slot = "data")
expr_after <- as.matrix(expr_after$SCT)
gsva.res_after <- gsva(expr_after, hallmark_genesets, method = "gsva")
#write.csv(gsva.res, "~/onedrive/Work/phD/phd_project/SiT/results/after_ln/after_ln_hallmark_gsva.csv")
rownames(gsva.res_after) <- gsub("HALLMARK_","",rownames(gsva.res))

colnames(gsva.res_after) <- c('Unknown','Fibroblast', 'Tcell Macrophage Fibroblast ',
                        'Epithelial cells // Unknown', 'Epithelial cells (CSTB+) and Macrophage',
                        'Malignant Epithelial cells , Macrophage, T cell','B cells, Tcells, Fibroblast,Macrophage',
                        'Epithelial cells (KRT17+,CSTB+) , Macrophage')
gsva.res_after <- gsva.res_after[,c('Unknown','Epithelial cells (KRT17+,CSTB+) , Macrophage',
                        'Epithelial cells // Unknown','Fibroblast','Tcell Macrophage Fibroblast ',
                        'B cells, Tcells, Fibroblast,Macrophage',
                        'Epithelial cells (CSTB+) and Macrophage','Malignant Epithelial cells , Macrophage, T cell'
                        
                        )]
### after_ln
expr_after_ln <- AverageExpression(after_ln, assays = "SCT", slot = "data")
expr_after_ln <- as.matrix(expr_after_ln$SCT)
gsva.res_after_ln <- gsva(expr_after_ln, hallmark_genesets, method = "gsva")
#write.csv(gsva.res, "~/onedrive/Work/phD/phd_project/SiT/results/after_ln/after_ln_hallmark_gsva.csv")
rownames(gsva.res_after_ln) <- gsub("HALLMARK_","",rownames(gsva.res))

colnames(gsva.res_after_ln) <- c("Tcells and Macrophage cluster 1","Fibroblast, Tcells and Macrophage",
                         "Tcells and Macrophage cluster 2","Tcells and Macrophage cluster 3",
                        "Epithelial cells and Tcells, Macrophage", "Malignant Epithelial cells, Tcells, Macrophage")
gsva.res_after_ln <- gsva.res_after_ln[,c("Tcells and Macrophage cluster 3","Epithelial cells and Tcells, Macrophage",
                        "Fibroblast, Tcells and Macrophage","Tcells and Macrophage cluster 1",
                        "Tcells and Macrophage cluster 2","Malignant Epithelial cells, Tcells, Macrophage")]


annotation_row <-  as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/hallmark.csv",col.names = F))
colnames(annotation_row) <- "Pathway"

pdf("~/onedrive/Work/phD/phd_project/SiT/results/before/before_hallmark.pdf",width = 6, height = 10)
pheatmap(gsva.res[rownames(annotation_row),],show_colnames = T, scale = "row", angle_col = "45", cluster_rows = F,cluster_cols = F,
                   color = viridis(50,option = "B"),annotation_row = annotation_row,gaps_row = c(6,17,24),
                   annotation_names_row = F)
dev.off()          

#### only Maliganant cells
gsva.res_malignant_cells <- cbind(gsva.res_before[,'Malignant Epithelial cells'],
                                  gsva.res_after[,'Malignant Epithelial cells , Macrophage, T cell'],
                                  gsva.res_after_ln[,"Malignant Epithelial cells, Tcells, Macrophage"])
colnames(gsva.res_malignant_cells) <- c('Before Malignant Epithelial cells', 'After Malignant Epithelial cells',
                                        'After LN Malignant Epithelial cells')

pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/scale_only_malignant_cells_hallmark.pdf",width = 10, height = 8)
pheatmap(gsva.res_malignant_cells[rownames(annotation_row),], show_colnames = T, angle_col = "45", cluster_rows = T,cluster_cols = F,
         col = viridis(10, option = "plasma" ),annotation_row = annotation_row,gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()  


## after down pathway
down_pathways <- c("HALLMARK_MITOTIC_SPINDLE","HALLMARK_HEDGEHOG_SIGNALING",
                   "HALLMARK_KRAS_SIGNALING_UP","HALLMARK_TGF_BETA_SIGNALING")
down_pathways_genes <- hallmark_genesets[down_pathways]
down_pathways_genes <- hallmark_genesets["HALLMARK_MITOTIC_SPINDLE"]
down_pathways_genes_all  <- unlist(down_pathways_genes)

# 提取原始表达矩阵
colnames(expr_before) <- c('Epithelial cells (KRT17+)','Fibroblast', 'Epithelial cells (CSTB+), T cells, B cells, Macrophage',
                                          'Bcells T cells','Malignant Epithelial cells', 'Epithelial cells (KRT17+, CSTB+), Macrophage')
expr_before_only_mag <- expr_before[, 'Malignant Epithelial cells']

colnames(expr_after) <- c('Unknown','Fibroblast', 'Tcell Macrophage Fibroblast ',
                        'Epithelial cells // Unknown', 'Epithelial cells (CSTB+) and Macrophage',
                        'Malignant Epithelial cells , Macrophage, T cell','B cells, Tcells, Fibroblast,Macrophage',
                        'Epithelial cells (KRT17+,CSTB+) , Macrophage')
expr_after_only_mag <- expr_after[, 'Malignant Epithelial cells , Macrophage, T cell']

colnames(expr_after_ln) <- c("Tcells and Macrophage cluster 3","Epithelial cells and Tcells, Macrophage",
                           "Fibroblast, Tcells and Macrophage","Tcells and Macrophage cluster 1",
                           "Tcells and Macrophage cluster 2","Malignant Epithelial cells, Tcells, Macrophage")
expr_after_ln_only_mag <- expr_after_ln[, "Malignant Epithelial cells, Tcells, Macrophage"]

total_genes <- intersect(intersect(names(expr_before_only_mag),names(expr_after_only_mag)),names(expr_after_ln_only_mag))
length(total_genes)
total_genes <- intersect(total_genes, down_pathways_genes_all)
length(total_genes)

combined_expr_only_mag <- cbind(expr_before_only_mag[total_genes],expr_after_only_mag[total_genes],
                                expr_after_ln_only_mag[total_genes])
colnames(combined_expr_only_mag) <- c('Before Malignant Epithelial cells', 'After Malignant Epithelial cells',
                                      'After LN Malignant Epithelial cells')

annotation_row <- data.frame(pathway = c(rep("MITOTIC_SPINDLE",length(down_pathways_genes[[1]]))),
                                         #rep("HEDGEHOG_SIGNALING",length(down_pathways_genes[[2]])),
                                         #rep("KRAS_SIGNALING_UP",length(down_pathways_genes[[3]])),
                                          #rep("TGF_BETA_SIGNALING",length(down_pathways_genes[[4]]))),
                             genes = c(down_pathways_genes[[1]]))
                                       #,down_pathways_genes[[2]],
                                       #down_pathways_genes[[3]],down_pathways_genes[[4]]))
annotation_row <- annotation_row[annotation_row$genes %in% total_genes,]
rownames(annotation_row) <- annotation_row$genes
pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/mitotic_spindle_pathway_genes_only_malignant_cells.pdf",width = 10, height = 50)
pheatmap(combined_expr_only_mag , show_colnames = T,show_rownames = T, 
         #scale = "row",
         angle_col = "45", cluster_rows = T,cluster_cols = F,
         col = viridis(10, option = "plasma" ),
         #annotation_row = annotation_row,
         #gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()  

#### Metadata ###3

  

