
suppressMessages(library(scRepertoire))
packageVersion("scRepertoire")
library(circlize)
library(scales)
####### 1. Data ########
sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/res_0.5.rds")
tcells_only_sce_anno_ori <- subset(sce_anno, subset = celltype %in% c("CD4+ Naive T", 
                                                                      "CD4+ Treg",
                                                                      "CD4+ Tfh",
                                                                      "CD8+ Tem",  
                                                                      "CD8+ Tcm", 
                                                                      "CD8+ Effctor T"))

####### 2. Combined with Tcell Seurat ######

before_file <- read.csv("~/onedrive/Work/phD/phd_project/SiT/rawdata/VDJ/before/filtered_contig_annotations.csv")
after_file <- read.csv("~/onedrive/Work/phD/phd_project/SiT/rawdata/VDJ/after/filtered_contig_annotations.csv")
after_ln_file <- read.csv("~/onedrive/Work/phD/phd_project/SiT/rawdata/VDJ/after_ln/filtered_contig_annotations.csv")

table(before_file$chain)
table(after_file$chain)
table(after_ln_file$chain)
contig_list <- list(before_file,after_file,after_ln_file)

combined <- combineTCR(contig_list, 
                       samples = c("before", "after", "after_ln"),
                       ID = c("before", "after", "after_ln"),
                       cells ="T-AB")
# 去除barcode前缀注释信息

combined$before_before$barcode <- gsub("before_before_","before_",combined$before_before$barcode)
combined$after_after$barcode <- gsub("after_after_","after_",combined$after_after$barcode)
combined$after_ln_after_ln$barcode <- gsub("after_ln_after_ln_","after_ln_",combined$after_ln_after_ln$barcode)


combined_tcr_sce <- combineExpression(combined, tcells_only_sce_anno_ori, 
                                      cloneCall="gene", 
                                      groupBy = "sample", 
                                      proportion = FALSE, 
                                      cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

####### 3. Visualization (Fig S13) #####
##### 3.1
pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_cloneType_dimplot.pdf",width = 15)
DimPlot(combined_tcr_sce, group.by = "cloneType") + scale_color_viridis(discrete = T)
DimPlot(combined_tcr_sce, group.by = "CTgene") + scale_color_viridis(discrete = T)

dev.off()
##### 3.2
# 分面展示不同样本的情况
pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_clonalOverlay.pdf",width = 30)

clonalOverlay(combined_tcr_sce, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, 
              facet = "orig.ident") + scale_color_viridis(discrete = T)

dev.off()
###### 3.3
pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_occupiedscRepertoire.pdf",width = 10)

occupiedscRepertoire(combined_tcr_sce, x.axis = "celltype") + scale_fill_viridis(discrete = T)
occupiedscRepertoire(combined_tcr_sce, x.axis = "orig.ident") + scale_fill_viridis(discrete = T)
dev.off()

###### 3.5


circles_1 <- getCirclize(combined_tcr_sce, groupBy = "celltype")
circles_2 <- getCirclize(combined_tcr_sce, groupBy = "orig.ident")

#Just assigning the normal colors to each cluster

#Graphing the chord diagram
mycolor <- viridis(6, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:6)]

pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_chordDiagram.pdf",width = 10)

circlize::chordDiagram(circles_1, self.link = 1, grid.col = mycolor) 
#circlize::chordDiagram(circles_2, self.link = 1, grid.col = grid.cols)+ scale_fill_viridis(discrete = T)
dev.off()

###### 3.6 
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))

StartracDiversity(combined_tcr_sce, type = "cloneType", 
                  sample = "orig.ident", by = "overall")

##### 3.7 
combined2 <- expression2List(combined_tcr_sce, group = "celltype")
combined3 <- expression2List(combined_tcr_sce, group = "orig.ident")

length(combined2) #now listed by cluster

pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_scores.pdf",width = 10,height = 5)

clonalDiversity(combined2, cloneCall = "nt") + scale_fill_viridis(discrete = T)
clonalDiversity(combined3, cloneCall = "nt") + scale_fill_viridis(discrete = T)

clonalHomeostasis(combined2, cloneCall = "nt") + scale_fill_viridis(discrete = T)
clonalHomeostasis(combined3, cloneCall = "nt") + scale_fill_viridis(discrete = T)

clonalProportion(combined2, cloneCall = "nt") + scale_fill_viridis(discrete = T)
clonalProportion(combined3, cloneCall = "nt") + scale_fill_viridis(discrete = T)

clonalOverlap(combined2, cloneCall="aa", method="overlap") + scale_color_viridis(discrete = T)
dev.off()
#### 3.8 
pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_compareClonotypes.pdf",width = 10,height = 5)

lengthContig(combined, cloneCall="aa", chains = "combined") + scale_fill_viridis(discrete = T)
lengthContig(combined, cloneCall="aa", chains = "TRA") + scale_fill_viridis(discrete = T)

compareClonotypes(combined, numbers = 10, 
                  samples = c("before_before", "after_after","after_ln_after_ln"), 
                  cloneCall="aa", graph = "alluvial") + scale_fill_viridis(discrete = T)
dev.off()

#### 3.9 
#设置group = "ID"，进行分组可视化
pdf("~/onedrive/Work/phD/phd_project/SiT/results/TCR/tcells_quantContig_abundanceContig.pdf",width = 10,height = 5)

quantContig(combined, cloneCall="gene", group = "ID", scale = TRUE)+ scale_fill_viridis(discrete = T)
abundanceContig(combined, cloneCall = "gene",group = "ID", scale = FALSE)+ scale_color_viridis(discrete = T)

abundanceContig(combined, cloneCall = "gene",group = "CTgene", scale = FALSE)+ scale_color_viridis(discrete = T)
abundanceContig(combined, cloneCall = "gene",group = "TCR1", scale = FALSE)+ scale_color_viridis(discrete = T)
dev.off()
