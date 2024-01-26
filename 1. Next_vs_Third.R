############# 0. Library import ########
library(cowplot)
library(ArchR)
############# 1. Percentage Illumina UMIs / Features seen in Nanopore data ################
###### data 
before@meta.data$nCount_profiled <- before@meta.data$nCount_ISOG*100/before@meta.data$nCount_Spatial
mean(before@meta.data$nCount_profiled)
# 34.42947
before@meta.data$nFeature_profiled <- before@meta.data

after@meta.data$nCount_profiled <- after@meta.data$nCount_ISOG*100/after@meta.data$nCount_Spatial
mean(after@meta.data$nCount_profiled)
# 25.14571
after@meta.data$nFeature_profiled <- after@meta.data

after_ln@meta.data$nCount_profiled <- after_ln@meta.data$nCount_ISOG*100/after_ln@meta.data$nCount_Spatial
mean(after_ln@meta.data$nCount_profiled)
# 22.99359
after_ln@meta.data$nFeature_profiled <- after_ln@meta.data

############# 2. Correlation illumina/Nanopore (r²) ##########
## before
p1 <- ggplot(data=before@meta.data, aes(x=nFeature_Spatial,y=nFeature_ISOG,fill = before@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("Features per cell") +
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore") +
  geom_rug(alpha = 0.08, position = "jitter")
x <- cor.test(before@meta.data$nFeature_Spatial, before@meta.data$nFeature_ISOG, method="pearson")$estimate
x*x
#0.9310187 
p2 <- ggplot(data=before@meta.data, aes(x=nCount_Spatial,y=nCount_ISOG,fill = before@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("UMIs per cell") + 
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore")+
  geom_rug(alpha = 0.08, position = "jitter")

x <- cor.test(before@meta.data$nCount_Spatial, before@meta.data$nCount_ISOG, method="pearson")$estimate
x*x
# 0.966456
## after 
p3 <- ggplot(data=after@meta.data, aes(x=nFeature_Spatial,y=nFeature_ISOG,fill = after@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("Features per cell") +
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore")+
  geom_rug(alpha = 0.08, position = "jitter")
x <- cor.test(after@meta.data$nFeature_Spatial, after@meta.data$nFeature_ISOG, method="pearson")$estimate
x*x
# 0.9527733
p4 <- ggplot(data=after@meta.data, aes(x=nCount_Spatial,y=nCount_ISOG,fill = after@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("UMIs per cell") + 
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore")+
  geom_rug(alpha = 0.08, position = "jitter")
x <- cor.test(after@meta.data$nCount_Spatial, after@meta.data$nCount_ISOG, method="pearson")$estimate
x*x
# 0.9900447
## after LN
p5 <- ggplot(data=after_ln@meta.data, aes(x=nFeature_Spatial,y=nFeature_ISOG,fill = after_ln@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("Features per cell") +
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore")+
  geom_rug(alpha = 0.08, position = "jitter")
x <- cor.test(after_ln@meta.data$nFeature_Spatial, after_ln@meta.data$nFeature_ISOG, method="pearson")$estimate
x*x
# 0.9628885
p6 <- ggplot(data=after_ln@meta.data, aes(x=nCount_Spatial,y=nCount_ISOG,fill = after_ln@active.ident)) +  geom_point(shape = 21, colour = "black") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle("UMIs per cell") + 
  #stat_cor(data=before@meta.data, aes(label = ..rr.label..),method = "pearson") +
  #geom_abline(slope=1, intercept=-15, linetype=3) +
  scale_fill_viridis(discrete = T) +
  labs(x="Illumina",y="Nanopore") +  ylim(0, 15000) +
  geom_rug(alpha = 0.08, position = "jitter") 

x <- cor.test(after_ln@meta.data$nCount_Spatial, after_ln@meta.data$nCount_ISOG, method="pearson")$estimate
x*x
# 0.9485484
pdf("results/illumina_vs_Nanopore_UMI_features.pdf",width = 15, height = 10)
ggarrange(p1, p3, p5,p2,p4,p6, ncol=3,nrow = 2,legend = "none")
dev.off()
############# 3. Correlation Illumina / Nanopore (r) ##########
### before
common.cells <- rownames(before@meta.data)
common.genes <- intersect(rownames(before@assays$Spatial@data), rownames(before@assays$ISOG@data))
common.genes.nomt <- common.genes[grep("mt-", common.genes, invert = TRUE)]
common.genes.nomtrp <- common.genes.nomt[grep("Rp[ls]", common.genes.nomt, invert = TRUE)]
length(common.genes.nomtrp)

df <- data.frame()
for (i in (1:length(common.cells))){
  c <- intersect(names(which(before@assays$Spatial@data[common.genes.nomtrp,common.cells[i]]>0)), names(which(before@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
  #print(paste(i,common.cells[i],length(c),sep=","))
  t <- data.frame(before@assays$Spatial@data[c,common.cells[i]], before@assays$ISOG@data[c,common.cells[i]])
  colnames(t) <- c("illu","nano")
  
  if(length(c)>100){
    df <- rbind(df, data.frame("Illu/ONT same cell",cor.test(t$illu, t$nano)$estimate))
  }
}
colnames(df) <- c("type","cor")
rownames(df) <- c()
median(df$cor, na.rm=TRUE)

dfmob.1 <- df

df2 <- data.frame()
for (i in (1:length(common.cells))){
  #print (i)
  random <- sample(1:length(common.cells), 2, replace=T)
  
  for (j in (1:length(random))){
    if(before@meta.data$celltype[i] != before@meta.data$celltype[random[j]]){
      c <- intersect(names(which(before@assays$Spatial@data[common.genes.nomtrp,common.cells[random[j]]]>0)), names(which(before@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
      t <- data.frame(before@assays$Spatial@data[c,common.cells[random[j]]], before@assays$ISOG@data[c,common.cells[i]])
      colnames(t) <- c("illu","nano")
      #print(paste(i,common.cells[i],length(c),sep=","))
      if(length(c)>100){
        df2 <- rbind(df2, data.frame("Illu/ONT different cell",cor.test(t$illu, t$nano)$estimate))
      }
    }
  }
}
colnames(df2) <- c("type","cor")
rownames(df2) <- c()
median(df2$cor, na.rm=TRUE)

dfmob.2 <- df2
# Nanopore/Illumina same cell
median(df$cor)
# Nanopore/Illumina different cell
median(df2$cor)
whole <- rbind(df,df2)
p1 <- ggplot(whole, aes(x=type, y=cor, fill=type)) +
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Transcriptomics correlation") +
  theme_classic()

### after
common.cells <- rownames(after@meta.data)
common.genes <- intersect(rownames(after@assays$Spatial@data), rownames(after@assays$ISOG@data))
common.genes.nomt <- common.genes[grep("mt-", common.genes, invert = TRUE)]
common.genes.nomtrp <- common.genes.nomt[grep("Rp[ls]", common.genes.nomt, invert = TRUE)]
length(common.genes.nomtrp)

df <- data.frame()
for (i in (1:length(common.cells))){
  c <- intersect(names(which(after@assays$Spatial@data[common.genes.nomtrp,common.cells[i]]>0)), names(which(after@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
  #print(paste(i,common.cells[i],length(c),sep=","))
  t <- data.frame(after@assays$Spatial@data[c,common.cells[i]], after@assays$ISOG@data[c,common.cells[i]])
  colnames(t) <- c("illu","nano")
  
  if(length(c)>100){
    df <- rbind(df, data.frame("Illu/ONT same cell",cor.test(t$illu, t$nano)$estimate))
  }
}
colnames(df) <- c("type","cor")
rownames(df) <- c()
median(df$cor, na.rm=TRUE)

dfmob.1 <- df

df2 <- data.frame()
for (i in (1:length(common.cells))){
  #print (i)
  random <- sample(1:length(common.cells), 2, replace=T)
  
  for (j in (1:length(random))){
    if(after@meta.data$celltype[i] != after@meta.data$celltype[random[j]]){
      c <- intersect(names(which(after@assays$Spatial@data[common.genes.nomtrp,common.cells[random[j]]]>0)), names(which(after@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
      t <- data.frame(after@assays$Spatial@data[c,common.cells[random[j]]], after@assays$ISOG@data[c,common.cells[i]])
      colnames(t) <- c("illu","nano")
      #print(paste(i,common.cells[i],length(c),sep=","))
      if(length(c)>100){
        df2 <- rbind(df2, data.frame("Illu/ONT different cell",cor.test(t$illu, t$nano)$estimate))
      }
    }
  }
}
colnames(df2) <- c("type","cor")
rownames(df2) <- c()
median(df2$cor, na.rm=TRUE)

dfmob.2 <- df2
# Nanopore/Illumina same cell
median(df$cor)
# Nanopore/Illumina different cell
median(df2$cor)
whole <- rbind(df,df2)
p2 <- ggplot(whole, aes(x=type, y=cor, fill=type)) +
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Transcriptomics correlation") +
  theme_classic()

### after_ln
common.cells <- rownames(after_ln@meta.data)
common.genes <- intersect(rownames(after_ln@assays$Spatial@data), rownames(after_ln@assays$ISOG@data))
common.genes.nomt <- common.genes[grep("mt-", common.genes, invert = TRUE)]
common.genes.nomtrp <- common.genes.nomt[grep("Rp[ls]", common.genes.nomt, invert = TRUE)]
length(common.genes.nomtrp)

df <- data.frame()
for (i in (1:length(common.cells))){
  c <- intersect(names(which(after_ln@assays$Spatial@data[common.genes.nomtrp,common.cells[i]]>0)), names(which(after_ln@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
  #print(paste(i,common.cells[i],length(c),sep=","))
  t <- data.frame(after_ln@assays$Spatial@data[c,common.cells[i]], after_ln@assays$ISOG@data[c,common.cells[i]])
  colnames(t) <- c("illu","nano")
  
  if(length(c)>100){
    df <- rbind(df, data.frame("Illu/ONT same cell",cor.test(t$illu, t$nano)$estimate))
  }
}
colnames(df) <- c("type","cor")
rownames(df) <- c()
median(df$cor, na.rm=TRUE)

dfmob.1 <- df

df2 <- data.frame()
for (i in (1:length(common.cells))){
  #print (i)
  random <- sample(1:length(common.cells), 2, replace=T)
  
  for (j in (1:length(random))){
    if(after_ln@meta.data$celltype[i] != after_ln@meta.data$celltype[random[j]]){
      c <- intersect(names(which(after_ln@assays$Spatial@data[common.genes.nomtrp,common.cells[random[j]]]>0)), names(which(after_ln@assays$ISOG@data[common.genes.nomtrp,common.cells[i]] > 0)))
      t <- data.frame(after_ln@assays$Spatial@data[c,common.cells[random[j]]], after_ln@assays$ISOG@data[c,common.cells[i]])
      colnames(t) <- c("illu","nano")
      #print(paste(i,common.cells[i],length(c),sep=","))
      if(length(c)>100){
        df2 <- rbind(df2, data.frame("Illu/ONT different cell",cor.test(t$illu, t$nano)$estimate))
      }
    }
  }
}
colnames(df2) <- c("type","cor")
rownames(df2) <- c()
median(df2$cor, na.rm=TRUE)

dfmob.2 <- df2
# Nanopore/Illumina same cell
median(df$cor)
# Nanopore/Illumina different cell
median(df2$cor)
whole <- rbind(df,df2)
p3 <- ggplot(whole, aes(x=type, y=cor, fill=type)) +
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Transcriptomics correlation") +
  theme_classic()

## combine
pdf("results/illumina_vs_Nanopore_transcriptomics_correlation.pdf", width=15, height=5, useDingbats=FALSE)
plot_grid(p1, p2, p3, ncol=3)
dev.off()


############# 4. ISO Distribution plots ########
## before
plot1 <- ggplot(data=before@meta.data, aes(nCount_ISO)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISO") + 
  theme_minimal() + 
  labs(x="nCount_ISO", y="Cells")

plot2 <- ggplot(data=before@meta.data, aes(nFeature_ISO)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISO") + 
  theme_minimal() + 
  labs(x="nFeature_ISO", y="Cells")

## after
plot3 <- ggplot(data=after@meta.data, aes(nCount_ISO)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISO") + 
  theme_minimal() + 
  labs(x="nCount_ISO", y="Cells")

plot4 <- ggplot(data=after@meta.data, aes(nFeature_ISO)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISO") + 
  theme_minimal() + 
  labs(x="nFeature_ISO", y="Cells")

## after_LN
plot5 <- ggplot(data=after_ln@meta.data, aes(nCount_ISO)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISO") + 
  theme_minimal() + 
  labs(x="nCount_ISO", y="Cells")

plot6 <- ggplot(data=after_ln@meta.data, aes(nFeature_ISO)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISO") + 
  theme_minimal() + 
  labs(x="nFeature_ISO", y="Cells")


pdf("results/illumina_vs_Nanopore_ISO_nCount_nFeature_distribution_plot.pdf",width = 15, height = 10)
plot_grid(plot1, plot3, plot5,plot2,plot4,plot6, ncol=3)
dev.off()

############# 5. ISO  celltype distribution plot #######
## nFeature_ISO
p1 <- ggplot(before@meta.data, aes(x = nFeature_ISO ,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C") +
  ggtitle( "Before")
p2 <- ggplot(after@meta.data, aes(x = nFeature_ISO ,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C")+
  ggtitle( "After")
p3 <- ggplot(after_ln@meta.data, aes(x = nFeature_ISO ,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C")+
  ggtitle( "After LN")
pdf("results/illumina_vs_Nanopore_ISO_nFeature_celltype_distribution_plot.pdf",width = 20, height = 5)
ggarrange(p1,p2,p3, ncol=3,legend ="none")
dev.off()

## nCounts ISO
p1 <- ggplot(before@meta.data, aes(x = nCount_ISO ,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C") +
  ggtitle( "Before")
p2 <- ggplot(after@meta.data, aes(x = nCount_ISO ,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C")+
  ggtitle( "After")
p3 <- ggplot(after_ln@meta.data, aes(x = nCount_ISO,y=celltype, fill=..y..))+
  geom_density_ridges_gradient(scale=3, rel_min_height=0.01, gradient_lwd = 1.)+
  scale_x_continuous(expand = c(0.01, 0))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="# of cells", option = "C")+
  ggtitle( "After LN")
pdf("results/illumina_vs_Nanopore_ISO_nCount_celltype_distribution_plot.pdf",width = 20, height = 5)
ggarrange(p1,p2,p3, ncol=3,legend ="none")
dev.off()


############# 6. ISOG Distribution plots ########
## before
plot1 <- ggplot(data=before@meta.data, aes(nCount_ISOG)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISOG") + 
  theme_minimal() + 
  labs(x="nCount_ISOG", y="Cells")

plot2 <- ggplot(data=before@meta.data, aes(nFeature_ISOG)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISOG") + 
  theme_minimal() + 
  labs(x="nFeature_ISOG", y="Cells")

## after
plot3 <- ggplot(data=after@meta.data, aes(nCount_ISOG)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISOG") + 
  theme_minimal() + 
  labs(x="nCount_ISOG", y="Cells")

plot4 <- ggplot(data=after@meta.data, aes(nFeature_ISOG)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISOG") + 
  theme_minimal() + 
  labs(x="nFeature_ISOG", y="Cells")

## after_LN
plot5 <- ggplot(data=after_ln@meta.data, aes(nCount_ISOG)) + 
  geom_histogram(breaks=seq(0, 16000, by = 250),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nCount_ISOG") + 
  theme_minimal() + 
  labs(x="nCount_ISOG", y="Cells")

plot6 <- ggplot(data=after_ln@meta.data, aes(nFeature_ISOG)) + 
  geom_histogram(breaks=seq(0, 4500, by = 125),col="black",aes(fill=nGene),fill="lightsteelblue3") + 
  scale_fill_gradient("Count", low = "green", high = "red") +
  labs(title="nFeature_ISOG") + 
  theme_minimal() + 
  labs(x="nFeature_ISOG", y="Cells")


pdf("results/illumina_vs_Nanopore_ISOG_nCount_nFeature_distribution_plot.pdf",width = 15, height = 10)
plot_grid(plot1, plot3, plot5,plot2,plot4,plot6, ncol=3)
dev.off()


############# 7. assignment efficiency ########
assignment_table <- read_excel("results/correction_table.xlsx")
assignment_table <- assignment_table[2:4,]
ordered_assignment_table <- data.frame(`assignment` = c(assignment_table$before,
                                                       assignment_table$after,assignment_table$after_LN),
                                       type = rep(assignment_table$Term,3),sample = c(rep("Before",3),rep("After",3),rep("After LN",3)))
assignment_eff <- sapply(strsplit(ordered_assignment_table$assignment,"\\("),"[[",2)
assignment_eff <- sapply(strsplit(assignment_eff,"\\%"),"[[",1)
ordered_assignment_table$assignment <- as.numeric(assignment_eff)
## plot
pdf("results/assignment_efficiency.pdf")
ggplot(ordered_assignment_table,aes(x= reorder(type,-assignment), y = assignment,fill = sample)) +
    geom_bar(stat="identity",position=position_dodge()) + scale_fill_brewer()  + xlab("")+
  ylab("% assignment")+
  geom_text(aes(label=as.character(assignment)),size=4,position=position_dodge(1))
dev.off()

############# 8. ISO UMAP vs SCT UMAP ############
## color palatte
my_color_palette <- hue_pal()(10)
names(my_color_palette) <- c(0,1:9)
####### before

DefaultAssay(before) <- "ISO"
before_iso_features <- FindVariableFeatures(before, selection.method = "vst", nfeatures = 10000)
before_iso_pca <- RunPCA(before_iso_features, features = VariableFeatures(object = before_iso_features))

use.pcs = 1:30
before_iso_cluster <- FindNeighbors(before_iso_pca, dims = use.pcs)
before_iso_cluster <- FindClusters(before_iso_cluster, resolution = c(1))
before_iso_cluster <- RunUMAP(before_iso_cluster, dims = use.pcs)
before_iso_cluster <- RunTSNE(before_iso_cluster, dims = use.pcs)

before_sct <- before
DefaultAssay(before_sct) <- "SCT"

pdf("results/before_iso_umap_tsne_compare.pdf")
DimPlot(before_iso_cluster, reduction = "umap", label = TRUE)
DimPlot(before_iso_cluster, reduction = "tsne", label = TRUE)
DimPlot(before_sct ,reduction = "umap", label = TRUE)
DimPlot(before_sct ,reduction = "tsne", label = TRUE)
dev.off()

# spatial cluster plot
before_sct <- FindClusters(before_sct, resolution = c(0.3))
before_iso_cluster <- FindClusters(before_iso_cluster, resolution = c(4))

## replace to same number
ori_idents <- before_iso_cluster$ISO_snn_res.2
loc_list <- list(which(ori_idents  == 0), which(ori_idents  == 1),which(ori_idents  == 2),which(ori_idents  == 3),
                 which(ori_idents  == 4),which(ori_idents  == 5),which(ori_idents  == 6),which(ori_idents  == 7),
                 which(ori_idents  == 8))
#rep_list <- c(6,4,2,3,0,5,1)
rep_list <- c(5,1,0,0,0,4,3,1,2)
for (i in 1:length(loc_list)){
  ori_idents[loc_list[[i]]] <- rep_list[i]
}
before_iso_cluster@active.ident <- ori_idents


pdf("results/before_iso_spatial_0.3_2.pdf")
SpatialDimPlot(before_sct,label = TRUE, cols =my_color_palette,label.size = 5,pt.size.factor = 3) 
SpatialDimPlot(before_iso_cluster,label = TRUE,cols =my_color_palette,  label.size = 5,pt.size.factor =3) 
dev.off()


####### after
DefaultAssay(after) <- "ISO"
after_iso_features <- FindVariableFeatures(after, selection.method = "vst", nfeatures = 10000)
after_iso_pca <- RunPCA(after_iso_features, features = VariableFeatures(object = after_iso_features))

use.pcs = 1:30
after_iso_cluster <- FindNeighbors(after_iso_pca, dims = use.pcs)
after_iso_cluster <- FindClusters(after_iso_cluster, resolution = c(1))
after_iso_cluster <- RunUMAP(after_iso_cluster, dims = use.pcs)
after_iso_cluster <- RunTSNE(after_iso_cluster, dims = use.pcs)

after_sct <- after
DefaultAssay(after_sct) <- "SCT"

pdf("results/after_iso_umap_tsne_compare.pdf")
DimPlot(after_iso_cluster, reduction = "umap", label = TRUE)
DimPlot(after_iso_cluster, reduction = "tsne", label = TRUE)
DimPlot(after_sct ,reduction = "umap", label = TRUE)
DimPlot(after_sct ,reduction = "tsne", label = TRUE)
dev.off()

##### spatial cluster plot
after_sct <- FindClusters(after_sct, resolution = c(0.5))
after_iso_cluster <- FindClusters(after_iso_cluster, resolution = c(0.9))

## replace to same number
ori_idents <- after_iso_cluster$ISO_snn_res.0.9
loc_list <- list(which(ori_idents  == 0), which(ori_idents  == 1),which(ori_idents  == 2),which(ori_idents  == 3),
              which(ori_idents  == 4),which(ori_idents  == 5),which(ori_idents  == 6),which(ori_idents  == 7),which(ori_idents  == 8))
#rep_list <- c(6,4,2,3,0,5,1)
rep_list <- c(4,1,3,5,0,1,6,2,0)
for (i in 1:length(loc_list)){
  ori_idents[loc_list[[i]]] <- rep_list[i]
}
after_iso_cluster@active.ident <- ori_idents

pdf("results/after_iso_spatial_0.5_0.9.pdf")
SpatialDimPlot(after_sct,label = TRUE, cols =my_color_palette,label.size = 5,pt.size.factor = 2) 
SpatialDimPlot(after_iso_cluster,label = TRUE,cols =my_color_palette,  label.size = 5,pt.size.factor = 2) 
dev.off()


###### after LN
DefaultAssay(after_ln) <- "ISO"
after_ln_iso_features <- FindVariableFeatures(after_ln, selection.method = "vst", nfeatures = 10000)
after_ln_iso_pca <- RunPCA(after_ln_iso_features, features = VariableFeatures(object = after_ln_iso_features))

use.pcs = 1:30
after_ln_iso_cluster <- FindNeighbors(after_ln_iso_pca, dims = use.pcs)
after_ln_iso_cluster <- FindClusters(after_ln_iso_cluster, resolution = c(1))
after_ln_iso_cluster <- RunUMAP(after_ln_iso_cluster, dims = use.pcs)
after_ln_iso_cluster <- RunTSNE(after_ln_iso_cluster, dims = use.pcs)

after_ln_sct <- after_ln
DefaultAssay(after_ln_sct) <- "SCT"

pdf("results/after_ln_iso_umap_tsne_compare.pdf")
DimPlot(after_ln_iso_cluster, reduction = "umap", label = TRUE)
DimPlot(after_ln_iso_cluster, reduction = "tsne", label = TRUE)
DimPlot(after_ln_sct ,reduction = "umap", label = TRUE)
DimPlot(after_ln_sct ,reduction = "tsne", label = TRUE)
dev.off()

##### spatial cluster plot
after_ln_sct <- FindClusters(after_ln_sct, resolution = c(0.5))
after_ln_iso_cluster <- FindClusters(after_ln_iso_cluster, resolution = c(0.7))

## replace to same number
ori_idents <- after_ln_iso_cluster$ISO_snn_res.0.7
loc_list <- list(which(ori_idents  == 0), which(ori_idents  == 1),which(ori_idents  == 2),which(ori_idents  == 3),
                 which(ori_idents  == 4),which(ori_idents  == 5))
rep_list <- c(0,3,1,2,5,4)
for (i in 1:length(loc_list)){
  ori_idents[loc_list[[i]]] <- rep_list[i]
}
after_ln_iso_cluster@active.ident <- ori_idents

pdf("results/after_ln_iso_spatial_0.5_0.7.pdf")
SpatialDimPlot(after_ln_sct,label = TRUE, cols =my_color_palette,label.size = 5,pt.size.factor = 1.3) 
SpatialDimPlot(after_ln_iso_cluster,label = TRUE,cols =my_color_palette,  label.size = 5,pt.size.factor = 1.3) 
dev.off()

##### save all isoform clustering results
saveRDS(after_ln_iso_cluster,"~/onedrive/Work/phD/phd_project/SiT/results/after_ln_isoform_clustering.rds")
saveRDS(after_iso_cluster,"~/onedrive/Work/phD/phd_project/SiT/results/after_isoform_clustering.rds")
saveRDS(before_iso_cluster,"~/onedrive/Work/phD/phd_project/SiT/results/before_isoform_clustering.rds")

############# 9.unintegrated illumina Nanopore spatial distribution #######
pdf("results/before_sct_spatial_four_method.pdf")
SpatialDimPlot(before_sct_1,label = TRUE, label.size = 5,pt.size.factor = 3)
SpatialDimPlot(before_sct_2, label = TRUE, label.size = 5,pt.size.factor = 3)
SpatialDimPlot(before_sct_3, label = TRUE, label.size = 5,pt.size.factor = 3)
dev.off()
## before
before_multi <- before
DefaultAssay(before_multi) <- "MULTI"
before_multi <- ScaleData(before_multi, features = rownames(sct_integrated))
before_multi_features <- FindVariableFeatures(before_multi, selection.method = "vst", nfeatures = 10000)
before_multi_pca <- RunPCA(before_multi_features, features = VariableFeatures(object = before_multi_features))

use.pcs = 1:30
before_multi_cluster <- FindNeighbors(before_multi_pca, dims = use.pcs)
before_multi_cluster <- FindClusters(before_multi_cluster, resolution = c(0.3))
before_multi_cluster <- RunUMAP(before_multi_cluster, dims = use.pcs)
before_multi_cluster <- RunTSNE(before_multi_cluster, dims = use.pcs)

before_sct <- before
DefaultAssay(before_sct) <- "SCT"
before_sct <- FindClusters(before_sct, resolution = c(0.1))
before_sct <- RunUMAP(before_sct, dims = use.pcs)
before_sct <- RunTSNE(before_sct, dims = use.pcs)

pdf("results/before_multi_umap_tsne_compare.pdf")
DimPlot(before_multi_cluster, reduction = "umap", label = TRUE)
DimPlot(before_multi_cluster, reduction = "tsne", label = TRUE)
DimPlot(before_sct ,reduction = "umap", label = TRUE)
DimPlot(before_sct ,reduction = "tsne", label = TRUE)
dev.off()

## spatial cluster plot
## 3 clusters
pdf("results/before_multi_spatial_3clusters.pdf")
SpatialDimPlot(before_sct,label = TRUE, label.size = 5,pt.size.factor = 3)
SpatialDimPlot(before_multi_cluster, label = TRUE, label.size = 5,pt.size.factor = 3)
dev.off()

before_sct <- FindClusters(before_sct, resolution = c(0.3))
before_multi_cluster <- FindClusters(before_multi_cluster, resolution = c(0.5))
pdf("results/before_multi_spatial_test.pdf")
SpatialDimPlot(before_sct,label = TRUE, label.size = 5,pt.size.factor = 3)
SpatialDimPlot(before_multi_cluster, label = TRUE, label.size = 5,pt.size.factor = 3)
dev.off()
### 效果一般 

## after 
after_multi <- after
DefaultAssay(after_multi) <- "MULTI"
after_multi <- ScaleData(after_multi, features = rownames(sct_integrated))
after_multi_features <- FindVariableFeatures(after_multi, selection.method = "vst", nfeatures = 2000)
after_multi_pca <- RunPCA(after_multi_features, features = VariableFeatures(object = after_multi_features))

use.pcs = 1:30
after_multi_cluster <- FindNeighbors(after_multi_pca, dims = use.pcs)
after_multi_cluster <- FindClusters(after_multi_cluster, resolution = c(1))
after_multi_cluster <- RunUMAP(after_multi_cluster, dims = use.pcs)
after_multi_cluster <- RunTSNE(after_multi_cluster, dims = use.pcs)

after_sct <- after
DefaultAssay(after_sct) <- "SCT"

pdf("results/after_multi_umap_tsne_compare.pdf")
DimPlot(after_multi_cluster, reduction = "umap", label = TRUE)
DimPlot(after_multi_cluster, reduction = "tsne", label = TRUE)
DimPlot(after_sct ,reduction = "umap", label = TRUE)
DimPlot(after_sct ,reduction = "tsne", label = TRUE)
dev.off()


