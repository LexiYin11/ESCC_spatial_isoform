######### 1.  提取表达矩阵 #######
expr_before <- AverageExpression(before, assays = "SCT", slot = "data")
expr_before <- as.matrix(expr_before$SCT)

xpr_after <- AverageExpression(after, assays = "SCT", slot = "data")
expr_after <- as.matrix(expr_after$SCT)

expr_after_ln <- AverageExpression(after_ln, assays = "SCT", slot = "data")
expr_after_ln <- as.matrix(expr_after_ln$SCT)

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

combined_expr_only_mag <- cbind(as.integer(expr_before_only_mag[total_genes]),as.integer(expr_after_only_mag[total_genes]),
                                as.integer(expr_after_ln_only_mag[total_genes]))
colnames(combined_expr_only_mag) <- c('Before Malignant Epithelial cells', 'After Malignant Epithelial cells',
                                      'After LN Malignant Epithelial cells')
rownames(combined_expr_only_mag)  <- total_genes

## coldata
coldata <- data.frame(ori_id = colnames(combined_expr_only_mag), group = c("Before", "After", "After_LN") )

plot.mat_diff = scale_mat(combined_expr_only_mag[rownames(dds.temp@assays@data$counts),],"row")
plot.mat_diff <- na.omit(plot.mat_diff )
dim(plot.mat_diff )
plot.mat_diff <- as.data.frame(plot.mat_diff)
######### 2.  selection #########
up_plot.mat_diff <- plot.mat_diff[0.5*plot.mat_diff$`After Malignant Epithelial cells` > plot.mat_diff$`Before Malignant Epithelial cells` & 0.5*plot.mat_diff$`After Malignant Epithelial cells`  > plot.mat_diff$`After LN Malignant Epithelial cells`,]
dim(up_plot.mat_diff)
down_plot.mat_diff <- plot.mat_diff[((0.05*plot.mat_diff$`After Malignant Epithelial cells`) <  plot.mat_diff$`Before Malignant Epithelial cells`) & ((0.05*plot.mat_diff$`After Malignant Epithelial cells`)  < plot.mat_diff$`After LN Malignant Epithelial cells`),]
dim(down_plot.mat_diff)

write.csv(up_plot.mat_diff,"~/onedrive/Work/phD/phd_project/SiT/results/function/up_malignant_genes.csv")
write.csv(down_plot.mat_diff,"~/onedrive/Work/phD/phd_project/SiT/results/function/down_malignant_genes.csv")

after_malignant_plot.mat_diff <- rbind(up_plot.mat_diff, down_plot.mat_diff)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/total_genes_only_malignant_cells.pdf",width = 10, height = 10)
pheatmap(after_malignant_plot.mat_diff  , show_colnames = T,show_rownames = F, 
         scale = "row",
         angle_col = "45", cluster_rows = T,cluster_cols = F,
         col = viridis(10, option = "plasma" ),
         #annotation_row = annotation_row,
         #gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()  

pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/down_genes_only_malignant_cells.pdf",width = 10, height = 20)
pheatmap(down_plot.mat_diff , show_colnames = T,show_rownames = T, 
         #scale = "row",
         angle_col = "45", cluster_rows = T,cluster_cols = F,
         col = viridis(10, option = "plasma" ),
         #annotation_row = annotation_row,
         #gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()  
##### 3. MSG score gene set only #####
msg_score_genes_df <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/function/final_msg_score_gene_df.csv")
msg_score_genes <- msg_score_genes_df$Approved.symbol

pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/msg_genes_malignant_cells_heatmap.pdf",width = 8, height = 8)
pheatmap(after_malignant_plot.mat_diff[msg_score_genes, ]  , show_colnames = T,show_rownames = T, 
         scale = "row",
         angle_col = "45", cluster_rows = T,cluster_cols = F,
         col = viridis(50, option = "plasma" ),
         #annotation_row = annotation_row,
         #gaps_row = c(6,17,24),
         annotation_names_row = F)

dev.off()  

##### GO ######

msg_scores_ego <- enrichGO(gene       = msg_score_genes_df$Ensembl.gene.ID ,
                OrgDb         = org.Hs.eg.db,
                keyType =       "ENSEMBL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
write.csv(msg_scores_ego ,"~/onedrive/Work/phD/phd_project/SiT/results/function/msg_genes_GO.csv")

msg_scores_ego@result <- msg_scores_ego@result[msg_scores_ego@result$Count > 4,]
msg_scores_ego@result <- msg_scores_ego@result[order(msg_scores_ego@result$p.adjust,decreasing = T),]
msg_scores_ego@result$p.adjust <- log10(msg_scores_ego@result$p.adjust )
dim(msg_scores_ego@result)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/msg_scores_GO_dotplot.pdf",w=5,h=5)

ggplot(msg_scores_ego@result , aes(x = p.adjust, y = reorder(Description, p.adjust))) + 
      geom_point(aes(size = Count,color = GeneRatio )) + theme_bw() + 
      scale_color_viridis(discrete = F,end = 0.9) +
      scale_size_continuous(range = c(3,6)) +
      theme(panel.grid = element_blank(),legend.position = "bottom")
dev.off()

allgene_pt <- pairwise_termsim(msg_scores_ego)
pdf("results/BaiduII_OV65/all_GO_FvsM_treeplot.pdf",w=10,h=10)
treeplot(allgene_pt)
dev.off()

pdf("results/BaiduII_OV65/all_GO_FvsM_emaplot.pdf",w=10,h=10)
emapplot(allgene_pt,cex_category=1.5,color = "p.adjust", layout = "kk")
dev.off()

#### KEGG #####
msg_gene_kegg <- enrichKEGG(gene       = msg_score_genes_df$NCBI.gene.ID ,
                      keyType =       "kegg",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 1)

dotplot(msg_gene_kegg )
####

