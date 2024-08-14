####  imm  vs others (Fig S8C)####
before$imm <- ifelse(before$seurat_clusters %in% c("2","3"),"Imm cells", "Others")
after$imm <- ifelse(after$seurat_clusters %in% c("2","5","6"),"Imm cells", "Others")
after_ln$imm <- "Imm cells"


Cellratio <- as.data.frame(rbind(prop.table(table(before$imm, before$orig.ident), margin = 2),
                   prop.table(table(after$imm, after$orig.ident), margin = 2),
                   prop.table(table(after_ln$imm, after_ln$orig.ident), margin = 1)))#计算各组样本不同细胞群比例
Cellratio$Sample <- c(rep(c("before","after"),each = 2), "after_ln")
Cellratio$Type <- c(rep(c("imm epithelial cells","Not imm epithelial cells"),2,each = 1),"imm epithelial cells")
colnames(Cellratio) <- c("Freq","Sample","Type")
library(ggplot2)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/ratio/imm_ratio.pdf")
ggplot(Cellratio) + 
  geom_bar(aes(x =Sample, y= Freq, fill = Type),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #scale_fill_viridis(option = "D",discrete = T)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()

####  malignant cells vs others (Fig S8B)####
before$malignant <- ifelse(before$celltype == 'malignant Epithelial cells',"malignant epithelial cells", "Not malignant epithelial cells")
after$malignant <- ifelse(after$celltype == 'malignant Epithelial cells , Macrophage, T cell', "malignant epithelial cells", "Not malignant epithelial cells")
after_ln$malignant <- ifelse(after_ln$celltype == "malignant Epithelial cells, Tcells, Macrophage", "malignant epithelial cells", "Not malignant epithelial cells")

as.data.frame(table(before$malignant))#查看各组细胞数
table(after$malignant)
table(after_ln$malignant)

Cellratio <- as.data.frame(rbind(prop.table(table(before$malignant, before$orig.ident), margin = 2),
                                 prop.table(table(after$malignant, after$orig.ident), margin = 2),
                                 prop.table(table(after_ln$malignant, after_ln$orig.ident), margin = 2)))#计算各组样本不同细胞群比例
Cellratio$sample <- rep(c("before","after", "after_ln"),each = 2)
Cellratio$Type <- rep(c("malignant epithelial cells","Not malignant epithelial cells"),3,each = 1)
colnames(Cellratio) <- c("Freq","Sample","Type")
library(ggplot2)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/ratio/malignant_ratio.pdf")
ggplot(Cellratio) + 
  geom_bar(aes(x =Sample, y= Freq, fill = Type),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #scale_fill_viridis(option = "D",discrete = T)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()
