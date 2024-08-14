####### 0. Library ###############
library(SPATA2)
library(Seurat)
library(viridis)
library(matrixStats)
library(ggstatsplot)
library(ggpubr)
library(dplyr)
library(STutility)
####### 1. Data import #######
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_ratio_labeled.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_ratio_labeled.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_ratio_labeled.rds")


tex_marker <- c("CD8A","CXCL13","CD3E","PDCD1")
tam_marker <- c("C1QC","LYZ","CD74","APOE","CSF1R")                

####### 2. Enrichment calculation (after_ln) ########
## Tex
after_ln_counts <- after_ln@assays$Spatial@data
after_ln_counts <- after_ln_counts[tex_marker,]
after_ln_counts <- as.data.frame(after_ln_counts)
after_ln_counts <- as.data.frame(t(after_ln_counts))
after_ln_tex_spots <- after_ln_counts %>% filter_all(all_vars(.>0) )

dim(after_ln_tex_spots)
dim(after_ln_counts)

## TAM
after_ln_counts <- after_ln@assays$Spatial@data
after_ln_counts <- after_ln_counts[tam_marker,]
after_ln_counts <- as.data.frame(after_ln_counts)
after_ln_counts <- as.data.frame(t(after_ln_counts))
after_ln_tam_spots <- after_ln_counts %>% filter_all(all_vars(.>0) )

dim(after_ln_tam_spots)
dim(after_ln_counts)

## label
after_ln$cell <- rownames(after_ln@meta.data)
after_ln$tex <- ifelse(after_ln$cell %in% rownames(after_ln_tex_spots), "CD8+CXCL13+Tex","No")
table(after_ln$tex)
after_ln$tam <- ifelse(after_ln$cell %in% rownames(after_ln_tam_spots), "C1QC+TAM","No")
table(after_ln$tam)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_colo.pdf")
SpatialDimPlot(after_ln, group.by = "tex")
SpatialDimPlot(after_ln, group.by = "tam")
SpatialDimPlot(after_ln, group.by = "colo")
dev.off()



####### 3. Enrichment calculation (after) ########
## Tex
after_counts <- after@assays$Spatial@data
after_counts <- after_counts[tex_marker,]
after_counts <- as.data.frame(after_counts)
after_counts <- as.data.frame(t(after_counts))
after_tex_spots <- after_counts %>% filter_all(all_vars(.>0) )

dim(after_tex_spots)
dim(after_counts)

## TAM
after_counts <- after@assays$Spatial@data
after_counts <- after_counts[tam_marker,]
after_counts <- as.data.frame(after_counts)
after_counts <- as.data.frame(t(after_counts))
after_tam_spots <- after_counts %>% filter_all(all_vars(.>0) )

dim(after_tam_spots)
dim(after_counts)

## label
after$cell <- rownames(after@meta.data)
after$tex <- ifelse(after$cell %in% rownames(after_tex_spots), "CD8+CXCL13+Tex","No")
table(after$tex)
after$tam <- ifelse(after$cell %in% rownames(after_tam_spots), "C1QC+TAM","No")
table(after$tam)
# after$colo <- ifelse(after$tex == "CD8+CXCL13+Tex" & 
#                           after$tam == "C1QC+TAM","Colocalization",
#                         ifelse(after$tex == "CD8+CXCL13+Tex" & 
#                                  after$tam != "C1QC+TAM","CD8+CXCL13+Tex",
#                                ifelse(after$tex != "CD8+CXCL13+Tex" & 
#                                         after$tam == "C1QC+TAM","C1QC+TAM","Others")))
# after$colo <- ifelse(after$tex == "CD8+CXCL13+Tex" & 
#                        after$tam == "C1QC+TAM","Colocalization","Others")
table(after$colo)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_colo.pdf")
SpatialDimPlot(after, group.by = "tex",pt.size.factor = 2.2)
SpatialDimPlot(after, group.by = "tam",pt.size.factor = 2.2)
SpatialDimPlot(after, group.by = "colo",pt.size.factor = 2.2)

dev.off()



####### 4. Enrichment calculation (before) ########
## Tex
before_counts <- before@assays$Spatial@data
before_counts <- before_counts[tex_marker,]
before_counts <- as.data.frame(before_counts)
before_counts <- as.data.frame(t(before_counts))
before_tex_spots <- before_counts %>% filter_all(all_vars(.>0) )

dim(before_tex_spots)
dim(before_counts)

## TAM
before_counts <- before@assays$Spatial@data
before_counts <- before_counts[tam_marker,]
before_counts <- as.data.frame(before_counts)
before_counts <- as.data.frame(t(before_counts))
before_tam_spots <- before_counts %>% filter_all(all_vars(.>0) )

dim(before_tam_spots)
dim(before_counts)

## label
before$cell <- rownames(before@meta.data)
before$tex <- ifelse(before$cell %in% rownames(before_tex_spots), "CD8+CXCL13+Tex","No")
table(before$tex)
before$tam <- ifelse(before$cell %in% rownames(before_tam_spots), "C1QC+TAM","No")
table(before$tam)
before$colo <- ifelse(before$tex == "CD8+CXCL13+Tex" & 
                       before$tam == "C1QC+TAM","Colocalization",
                     ifelse(before$tex == "CD8+CXCL13+Tex" & 
                              before$tam != "C1QC+TAM","CD8+CXCL13+Tex",
                            ifelse(before$tex != "CD8+CXCL13+Tex" & 
                                     before$tam == "C1QC+TAM","C1QC+TAM","Others")))
before$colo <- ifelse(before$tex == "CD8+CXCL13+Tex" & 
                       before$tam == "C1QC+TAM","Colocalization","Others")
table(before$colo)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/before_colo.pdf")
SpatialDimPlot(before, group.by = "tex",pt.size.factor = 2.2)
SpatialDimPlot(before, group.by = "tam",pt.size.factor = 2.2)
SpatialDimPlot(before, group.by = "colo",pt.size.factor = 2.2)

dev.off()

####### 5. Ratio label (Fig 4G,N) #####
## before
before_filter <- subset(before, subset = ratio !=  "Inf")
#before_filter <- subset(before, subset = is.nan(ratio1))
before_filter$ratio = before_filter@assays$RATIO@data["ratio1",]
before_filter_metadata <- before_filter@meta.data
before_filter_metadata <- before_filter_metadata[!is.nan(before_filter_metadata$ratio),]
dim(before_filter_metadata)

## after_ln
## filter all NaN Inf
after_ln_filter <- subset(after_ln, subset = RATIO !=  "Inf")
#after_ln_filter <- subset(after_ln, subset = is.nan(ratio1))
after_ln_filter$ratio = after_ln_filter@assays$RATIO@data["ratio1",]
after_ln_filter_metadata <- after_ln_filter@meta.data
after_ln_filter_metadata <- after_ln_filter_metadata[!is.nan(after_ln_filter_metadata$ratio),]
dim(after_ln_filter_metadata)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_ratio_total_enrichment.pdf")

ggboxplot( after_ln_filter_metadata,x = "total_enrichment", y = "ratio" ,fill = "total_enrichment") +
  stat_compare_means() + ylab("CD74-202 / CD74-201") +ggtitle("AfterLN")
dev.off()

## after
## filter all NaN Inf
after_filter <- subset(after, subset = ratio !=  "Inf")
#after_filter <- subset(after, subset = is.nan(ratio1))
after_filter$ratio = after_filter@assays$RATIO@data["ratio1",]
after_filter_metadata <- after_filter@meta.data
after_filter_metadata <- after_filter_metadata[!is.nan(after_filter_metadata$ratio),]
dim(after_filter_metadata)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ratio_total_enrichment.pdf")
ggboxplot( after_filter_metadata,x = "total_enrichment", y = "ratio" ,fill = "total_enrichment") +
  stat_compare_means() + ylab("CD74-202 / CD74-201") +ggtitle("After")
dev.off()



####### 6. Find Neighbor ######

## Function
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
cell <- rownames(after@meta.data)[155]
coordinates <- as.data.frame(coordinates)
cell_neighbors <- get_neighbor_function(as.data.frame(coordinates), cell)

after$label_test <- ifelse(after$cell %in% cell_neighbors, "neighbors",
                           ifelse(after$cell == cell, "cell", "others"))
table(after$label_test)                              
SpatialDimPlot(after, group.by = "label_test")

saveRDS(after, "~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ratio_labeled.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ratio_labeled.rds")
saveRDS(after_ln, "~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_ratio_labeled.rds")



####### 7. label neighbor enrichment ####
## Functions
# label method
enrichment_function_list <- function(seu_obj, cell_list, type,coordinates){
  enrichment_label_list <- c()
  metadata <- seu_obj@meta.data
  for (cell in cell_list){
    label <- metadata[cell,type]
    enrichment_label_list <- c(enrichment_label_list, 
                               enrichment_function(seu_obj, cell, label, coordinates))
  }
  return(enrichment_label_list)
}
enrichment_function <- function(seu_obj, cell, type, coordinates){
  metadata <- seu_obj@meta.data
  if (type == "CD8+CXCL13+Tex"){
    cell_neighbors<- get_neighbor_function(coordinates,cell=cell)
    all_cells <- c(cell_neighbors, cell)
    cell_neighbors_meta <- metadata[metadata$cell %in% all_cells,]
    if ("C1QC+TAM" %in% cell_neighbors_meta$tam){
      return("co-enrichment")
    } else {
      return("No")
    }
  } else if (type == "C1QC+TAM"){
    cell_neighbors<- get_neighbor_function(coordinates,cell=cell)
    all_cells <- c(cell_neighbors, cell)
    cell_neighbors_meta <- metadata[metadata$cell %in% all_cells,]
    if ("CD8+CXCL13+Tex" %in% cell_neighbors_meta$tex){
      return("co-enrichment")
    } else {
      return("No")
    }
  } else if (type == "No"){
    return("No")
  }
  
}

# score method
enrichment_function_list_score <- function(seu_obj, cell_list, type, coordinates){
  enrichment_label_list <- c()
  metadata <- seu_obj@meta.data
  for (cell in cell_list){
    label <- metadata[cell,type]
    enrichment_label_list <- c(enrichment_label_list, 
                               enrichment_function_score(seu_obj, cell, label, coordinates))
  }
  return(enrichment_label_list)
}
enrichment_function_score <- function(seu_obj, cell, type, coordinates){
  metadata <- seu_obj@meta.data
  if (type == "CD8+CXCL13+Tex"){
    cell_neighbors<- get_neighbor_function(coordinates,cell=cell)
    all_cells <- c(cell_neighbors, cell)
    cell_neighbors_meta <- metadata[metadata$cell %in% all_cells,]
    if ("C1QC+TAM" %in% cell_neighbors_meta$tam){
      return(sum(cell_neighbors_meta$tam == "C1QC+TAM"))
    } else {
      return(0)
    }
  } else if (type == "C1QC+TAM"){
    cell_neighbors<- get_neighbor_function(coordinates,cell=cell)
    all_cells <- c(cell_neighbors, cell)
    cell_neighbors_meta <- metadata[metadata$cell %in% all_cells,]
    if ("CD8+CXCL13+Tex" %in% cell_neighbors_meta$tex){
      return(sum(cell_neighbors_meta$tex == "CD8+CXCL13+Tex"))
    } else {
      return(0)
    }
  } else if (type == "No"){
    return(0)
  }
  
}

##### 7.1 after_ln
after_ln_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after_LN/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
## label tex enrichment 
after_ln$tex_enrichment <- enrichment_function_list(after_ln,after_ln$cell, "tex",after_ln_coordinates )
table(after_ln$tex_enrichment)


## label tam enrichment 
after_ln$tam_enrichment <- enrichment_function_list(after_ln,after_ln$cell, "tam",after_ln_coordinates )
table(after_ln$tam_enrichment)

## combine label
after_ln$total_enrichment <- ifelse(after_ln$tex_enrichment == "co-enrichment" | after_ln$tam_enrichment == "co-enrichment",
                                    "co-enrichment", "Others")
table(after_ln$total_enrichment)


## score tex enrichment 
after_ln$tex_enrichment_score <- enrichment_function_list_score(after_ln,after_ln$cell, "tex", after_ln_coordinates)
table(after_ln$tex_enrichment_score)


## score tam enrichment 
after_ln$tam_enrichment_score <- enrichment_function_list_score(after_ln,after_ln$cell, "tam" ,after_ln_coordinates)
table(after_ln$tam_enrichment_score) 

## combine score
# keep the highest enrichment score or add up
after_ln$total_enrichment_score <- ifelse(after_ln$tex_enrichment_score == 0, after_ln$tam_enrichment_score,
                                          ifelse(after_ln$tam_enrichment_score == 0, after_ln$tex_enrichment_score,
                                                 after_ln$tam_enrichment_score + after_ln$tex_enrichment_score ))
table(after_ln$total_enrichment_score)
#table(after_ln$colo)
SpatialFeaturePlot(after_ln, features = "total_enrichment_score")

####### 7.2 after
after_coordinates <-  LoadSpatialCoordinates("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_After/filtered_feature_bc_matrix/spatial/tissue_positions_list.csv") |>
  select(all_of(c("barcode", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres", "sampleID")))
## label tex enrichment 
after$tex_enrichment <- enrichment_function_list(after,after$cell, "tex",after_coordinates )
table(after$tex_enrichment)


## label tam enrichment 
after$tam_enrichment <- enrichment_function_list(after,after$cell, "tam",after_coordinates )
table(after$tam_enrichment)

## combine label
after$total_enrichment <- ifelse(after$tex_enrichment == "co-enrichment" | after$tam_enrichment == "co-enrichment",
                                    "co-enrichment", "Others")
table(after$total_enrichment)


## score tex enrichment 
after$tex_enrichment_score <- enrichment_function_list_score(after,after$cell, "tex", after_coordinates)
table(after$tex_enrichment_score)
table(after$tex_ratio_label)
#table(after$colo)

## score tam enrichment 
after$tam_enrichment_score <- enrichment_function_list_score(after,after$cell, "tam" ,after_coordinates)
table(after$tam_enrichment_score) 

## combine score
# keep the highest enrichment score or add up
after$total_enrichment_score <- ifelse(after$tex_enrichment_score == 0, after$tam_enrichment_score,
                                          ifelse(after$tam_enrichment_score == 0, after$tex_enrichment_score,
                                                 after$tam_enrichment_score + after$tex_enrichment_score ))
table(after$total_enrichment_score)


SpatialFeaturePlot(after, features = "total_enrichment_score")

####### 8. Ratio to enrichment label (Fig 4E,L) #####
## 8.1 after_ln
after_ln_filter <- subset(after_ln, subset = ratio !=  "Inf")
summary(after_ln_filter$ratio)
after_ln_filter$ratio_label <-   ifelse(after_ln_filter$ratio <= 1, "low",
                                      ifelse(after_ln_filter$ratio > 1 & after_ln_filter$ratio <= 15 , "medium",
                                             ifelse(after_ln_filter$ratio > 15, "high", NA)))
after_ln_filter_metadata <- after_ln_filter@meta.data
after_ln_filter_metadata <- after_ln_filter_metadata[!is.nan(after_ln_filter_metadata$ratio),]
after_ln_filter_metadata$ratio_label <- factor(after_ln_filter_metadata$ratio_label, levels = c("low", "medium", "high"))
table(after_ln_filter_metadata$ratio_label)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_ratio_to_enrichment_points.pdf")

ggbarstats(after_ln_filter_metadata, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 
after_ln_filter_metadata_1 <- after_ln_filter_metadata[after_ln_filter_metadata$ratio_label %in% c("low","medium"),]
ggbarstats(after_ln_filter_metadata_1, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

after_ln_filter_metadata_2 <- after_ln_filter_metadata[after_ln_filter_metadata$ratio_label %in% c("low","high"),]
ggbarstats(after_ln_filter_metadata_2, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

after_ln_filter_metadata_3 <- after_ln_filter_metadata[after_ln_filter_metadata$ratio_label %in% c("high","medium"),]
ggbarstats(after_ln_filter_metadata_3, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

dev.off()

## 8.2 after
after_filter <- subset(after, subset = ratio !=  "Inf")
summary(after_filter$ratio)
after_filter$ratio_label <- ifelse(after_filter$ratio <= 1, "low",
                                      ifelse(after_filter$ratio > 1 & after_filter$ratio <= 15 , "medium",
                                             ifelse(after_filter$ratio > 15, "high", NA)))
after_filter_metadata <- after_filter@meta.data
after_filter_metadata <- after_filter_metadata[!is.nan(after_filter_metadata$ratio),]
after_filter_metadata$ratio_label <- factor(after_filter_metadata$ratio_label, levels = c("low", "medium", "high"))
table(after_filter_metadata$ratio_label)
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ratio_to_enrichment.pdf")

ggbarstats(after_filter_metadata, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 
after_filter_metadata_1 <- after_filter_metadata[after_filter_metadata$ratio_label %in% c("low","medium"),]
ggbarstats(after_filter_metadata_1, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

after_filter_metadata_2 <- after_filter_metadata[after_filter_metadata$ratio_label %in% c("low","high"),]
ggbarstats(after_filter_metadata_2, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

after_filter_metadata_3 <- after_filter_metadata[after_filter_metadata$ratio_label %in% c("high","medium"),]
ggbarstats(after_filter_metadata_3, total_enrichment, ratio_label, ggstatsplot.layer = FALSE, title = "AfterLN") 

dev.off()

### 8.3 Pie chart 
## function
piechart_group <- function(Cellratio, group){
  cur_df <- Cellratio[Cellratio$Var2 == group,]
  cur_df$mylabel <- paste0(round(cur_df$Freq,4) * 100,"%")
  pie = ggplot(cur_df, aes(x="", y=Freq, fill=Var1)) + geom_bar(stat="identity", width=1)+
    scale_fill_manual(values = c("#E76D29", "#370857"))
  
  
  # Convert to pie (polar coordinates) and add labels
  pie = pie + coord_polar("y", start=0) + geom_text(aes(label = mylabel), position = position_stack(vjust = 0.5))
  
  # Remove labels and add title
  pie = pie + labs(x = NULL, y = NULL, fill = NULL)
  
  # Tidy up the theme
  pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank(),
                                      plot.title = element_text(hjust = 0.5, color = "#666666"))
  return(pie)
}

## after
ratio_after <- as.data.frame(prop.table(table(after_filter_metadata$total_enrichment,after_filter_metadata$ratio_label ), margin = 2))#计算各组样本不同细胞群比例

p1 <- piechart_group(ratio_after, "low")
p2 <- piechart_group(ratio_after, "medium")
p3 <-  piechart_group(ratio_after, "high")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_enrichment_pie.pdf")
ggarrange(p1, p2,p3, common.legend = T, ncol = 3)
dev.off()

## after_ln
ratio_after_ln <- as.data.frame(prop.table(table(after_ln_filter_metadata$total_enrichment,after_ln_filter_metadata$ratio_label ), margin = 2))#计算各组样本不同细胞群比例

p1 <- piechart_group(ratio_after_ln, "low")
p2 <- piechart_group(ratio_after_ln, "medium")
p3 <-  piechart_group(ratio_after_ln, "high")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_enrichment_pie.pdf")
ggarrange(p1, p2,p3, common.legend = T, ncol = 3)
dev.off()
####### 9. Ratio to enrichment score (Fig 4F,M) ######
## after_ln
after_ln_filter <- subset(after_ln, subset = ratio !=  "Inf")
summary(after_ln_filter$ratio)
after_ln_filter$ratio_label <- ifelse(after_ln_filter$ratio <= 1, "low",
                                      ifelse(after_ln_filter$ratio > 1 & after_ln_filter$ratio <= 15 , "medium",
                                             ifelse(after_ln_filter$ratio > 15, "high", NA)))
after_ln_filter_metadata <- after_ln_filter@meta.data
after_ln_filter_metadata <- after_ln_filter_metadata[!is.nan(after_ln_filter_metadata$ratio),]
after_ln_filter_metadata$ratio_label <- factor(after_ln_filter_metadata$ratio_label, levels = c("low","medium", "high"))


pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_ratio_to_enrichment_boxplot.pdf")
ggboxplot(after_ln_filter_metadata,x = "ratio_label", y = "total_enrichment_score" ,fill = "ratio_label") + 
  stat_compare_means(comparisons = list(c("low", "medium"), c("medium", "high"))) + scale_fill_manual(values = c("#1B042D","#992963","#E76D29" ))
dev.off()

## after
after_filter <- subset(after, subset = ratio !=  "Inf")
summary(after_filter$ratio)
after_filter$ratio_label <- ifelse(after_filter$ratio <= 1, "low",
                                      ifelse(after_filter$ratio > 1 & after_filter$ratio <= 15 , "medium",
                                             ifelse(after_filter$ratio > 15, "high", NA)))
after_filter_metadata <- after_filter@meta.data
after_filter_metadata <- after_filter_metadata[!is.nan(after_filter_metadata$ratio),]
after_filter_metadata$ratio_label <- factor(after_filter_metadata$ratio_label, levels = c("low","medium", "high"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ratio_to_enrichment_boxplot.pdf")
ggboxplot(after_filter_metadata,x = "ratio_label", y = "total_enrichment_score" ,fill = "ratio_label") + 
  stat_compare_means(comparisons = list(c("low", "medium"), c("medium", "high"))) + scale_fill_manual(values = c("#1B042D","#992963","#E76D29" ))
dev.off()

####### 10. Visualization ratio and co-enrichment (Fig 4D,K) ####
treatment_col = viridis(6,option = "B",direction = -1)
treatment_col <- colorRampPalette(treatment_col)(13)
SpatialColors <- rev(x = brewer.pal(n = 11, name = "Spectral"))
scales::show_col(treatment_col)

before$ratio_label <- ifelse(before$ratio == Inf, "Others",
                            ifelse(before$ratio <= quantile(counts, 0.25), "Low",
                                   ifelse(before$ratio > quantile(counts, 0.25) & before$ratio <= quantile(counts, 0.75) , "Medium",
                                          ifelse(before$ratio > quantile(counts, 0.75), "High", "Others"))))
before$ratio_label <- ifelse(is.na(before$ratio_label), "Others",before$ratio_label)
table(before$ratio_label)
after$ratio_label <- ifelse(after$ratio == Inf, "Others",
                            ifelse(after$ratio <= quantile(counts, 0.25), "Low",
                                   ifelse(after$ratio > quantile(counts, 0.25) & after$ratio <= quantile(counts, 0.75) , "Medium",
                                          ifelse(after$ratio > quantile(counts, 0.75), "High", "Others"))))
after$ratio_label <- ifelse(is.na(after$ratio_label), "Others",after$ratio_label)

after_ln$ratio_label <- ifelse(after_ln$ratio == Inf, "Others",
                               
                               ifelse(after_ln$ratio <= quantile(counts, 0.25), "Low",
                            ifelse(after_ln$ratio > quantile(counts, 0.25) & after_ln$ratio <= quantile(counts, 0.75) , "Medium",
                                   ifelse(after_ln$ratio > quantile(counts, 0.75), "High", "Others"))))
after_ln$ratio_label <- ifelse(is.na(after_ln$ratio_label), "Others",after_ln$ratio_label)

before$ratio_label <- factor(before$ratio_label, levels = c("High", "Medium","Low","Others"))

after$ratio_label <- factor(after$ratio_label, levels = c("High", "Medium","Low","Others"))
after_ln$ratio_label <- factor(after_ln$ratio_label, levels = c("High", "Medium","Low","Others"))

pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/before_enrichment.pdf")
SpatialDimPlot(before, group.by = "tex", pt.size.factor = 3)+ scale_fill_manual(values = c("#F49015", "#1B042D"))
SpatialDimPlot(before, group.by = "tam", pt.size.factor = 3)+ scale_fill_manual(values = c("#F49015", "#1B042D"))
#SpatialDimPlot(after_ln, group.by = "total_enrichment", pt.size.factor = 3) + scale_fill_manual(values = c("#E76D29", "#1B042D"))
SpatialDimPlot(before, group.by = "ratio_label", pt.size.factor = 3)  + scale_fill_manual(values = c("#FCD963","#E76D29","#992963", "#1B042D"))
SpatialFeaturePlot(before, features="ratio", pt.size.factor = 3) + scale_fill_gradientn(colors= SpatialColors,limits = c(0, 10))
#SpatialFeaturePlot(after_ln, features="total_enrichment_score", pt.size.factor = 1.5)

dev.off()


pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_ln_enrichment.pdf")
SpatialDimPlot(after_ln, group.by = "tex", pt.size.factor = 1.5)+ scale_fill_manual(values = c("#F49015", "#1B042D"))
SpatialDimPlot(after_ln, group.by = "tam", pt.size.factor = 1.5)+ scale_fill_manual(values = c("#F49015", "#1B042D"))
SpatialDimPlot(after_ln, group.by = "total_enrichment", pt.size.factor = 1.5) + scale_fill_manual(values = c("#E76D29", "#1B042D"))
SpatialDimPlot(after_ln, group.by = "ratio_label", pt.size.factor = 1.5)  + scale_fill_manual(values = c("#FCD963","#E76D29","#992963", "#1B042D"))
SpatialFeaturePlot(after_ln, features="ratio", pt.size.factor = 1.5) + scale_fill_gradientn(colors= SpatialColors,limits = c(0, 10))
SpatialFeaturePlot(after_ln, features="total_enrichment_score", pt.size.factor = 1.5)

dev.off()

pdf("~/onedrive/Work/phD/phd_project/SiT/results/shijie/after_enrichment.pdf")
SpatialDimPlot(after, group.by = "tex", pt.size.factor = 2.2) + scale_fill_manual(values = c("#F49015", "#370857"))
SpatialDimPlot(after, group.by = "tam", pt.size.factor = 2.2) + scale_fill_manual(values = c("#F49015", "#370857"))
SpatialDimPlot(after, group.by = "total_enrichment", pt.size.factor = 2.2) + scale_fill_manual(values = c("#E76D29", "#370857"))
SpatialDimPlot(after, group.by = "ratio_label", pt.size.factor = 2.2) + scale_fill_manual(values = c("#FCD963","#E76D29","#992963", "#1B042D"))
SpatialFeaturePlot(after, features="ratio", pt.size.factor = 2.2) + scale_fill_gradientn(colors= SpatialColors,limits = c(0, 30))
SpatialFeaturePlot(after, features="total_enrichment_score", pt.size.factor = 2.2)

dev.off()

table(after$ratio_label)
table(after_ln$ratio_label)
SpatialFeaturePlot(after, features="ratio", pt.size.factor = 2.2)
SpatialFeaturePlot(after, features="total_enrichment_score", pt.size.factor = 2.2)

