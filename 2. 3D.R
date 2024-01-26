############# 0. library import ##############
library("imager")
library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(biomaRt)
library(STutility)
############# 1. data prep #######
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_added_isoform.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_added_isoform.rds")
after_ln <- readRDS("/Users/lexiiiii/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_added_isoform.rds")

### three samples 3d data prepration 
list_3d <- c()
file_list <- c("~/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/1.Count/P1_before/filtered_feature_bc_matrix",
               "~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after/filtered_feature_bc_matrix",
               "~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002/1.Count/p1_after_ln/filtered_feature_bc_matrix")

seurat_list <- c(before,after,after_ln)
for (i in 1:length(seurat_list)){
  file <- file_list[i]
  seurat <- seurat_list[[i]]
  setwd(file)
  ###### 3 slices from the beginning 
  samples <- c("filtered_feature_bc_matrix.h5","filtered_feature_bc_matrix.h5","filtered_feature_bc_matrix.h5")
  imgs <- c("spatial/tissue_hires_image.png","spatial/tissue_hires_image.png","spatial/tissue_hires_image.png")
  spotfiles <- c("spatial/tissue_positions_list.csv","spatial/tissue_positions_list.csv","spatial/tissue_positions_list.csv")
  json <- c("spatial/scalefactors_json.json","spatial/scalefactors_json.json","spatial/scalefactors_json.json")
  infoTable <- data.frame(samples, imgs, spotfiles, json, stringsAsFactors = FALSE)
  
  se <- InputFromTable(infotable = infoTable, 
                       min.gene.count = 0, 
                       min.gene.spots = 0,
                       min.spot.count = 0,
                       platform =  "Visium")
  
  se <- LoadImages(se, time.resolve = F, verbose = T)
  se <- MaskImages(object = se, compactness = 2)
  ImagePlot(se, ncols = 2, method = "raster", type = "masked", darken = F) # Masked image
  se <- AlignImages(se)

  ##### Add iso and change isoform level into same name 
  ## RPS9 have three isoform so change their names into RPS9 and delete other two isoforms 
  isoforms_RPS9 <- c("RPS9..ENST00000391751","RPS9..ENST00000391753","RPS9..ENST00000302907")
  
  df1 <- as.data.frame(t(seurat@assays$ISO@counts))
  #df1 <- df1[,isoforms_RPS9[1]]
  colnames(df1)[which(colnames(df1) == isoforms_RPS9[1])] <- "RPS9"
  df1 <- df1[,!colnames(df1) %in% isoforms_RPS9[2:3]]
  rownames(df1) <- paste0(rownames(seurat@meta.data),"_1",sep="")
  dim(df1)
  
  df2 <- as.data.frame(t(seurat@assays$ISO@counts))
  colnames(df2)[which(colnames(df2) == isoforms_RPS9[2])] <- "RPS9"
  df2 <- df2[,!colnames(df2) %in% isoforms_RPS9[c(1,3)]]
  rownames(df2) <- paste0(rownames(seurat@meta.data),"_2",sep="")
  dim(df2)
  
  df3 <- as.data.frame(t(seurat@assays$ISO@counts))
  colnames(df3)[which(colnames(df3) == isoforms_RPS9[3])] <- "RPS9"
  df3 <- df3[,!colnames(df3) %in% isoforms_RPS9[c(1,2)]]
  rownames(df3) <- paste0(rownames(seurat@meta.data),"_3",sep="")
  dim(df3)
  
  dd <- rbind(df1,df2,df3)
  new_cells <- setdiff(rownames(se@meta.data),rownames(dd))
  new_df <- dd[1:length(new_cells),]
  rownames(new_df) <- new_cells
  dim(new_df)
  dd <- rbind(dd,new_df)
  se[["ISO"]] <- CreateAssayObject(data = t(dd))
  se <- NormalizeData(object = se, assay = "ISO")
  se <- ScaleData(se, assay = "ISO")
  se <- MaskImages(object = se)
  DefaultAssay(se) <- "ISO"
  se_3d <- Create3DStack(se)
  stack_3d <- setNames(GetStaffli(se_3d)@scatter.data, c("x", "y", "z", "grid.cell"))
  list_3d <- c(list_3d,se_3d)
}

############# 2. Plotting #######
## Feature plot
setwd("~/onedrive/Work/phD/phd_project/SiT")

interpolated.data.1 <- FeaturePlot3D(list_3d[[1]], features = "RPS9", return.data = TRUE)
interpolated.data.2 <- FeaturePlot3D(list_3d[[2]], features = "RPS9", return.data = TRUE)
interpolated.data.3 <- FeaturePlot3D(list_3d[[3]], features = "RPS9", return.data = TRUE)

df <- data.frame(value = c(interpolated.data.1$val,interpolated.data.2$val,interpolated.data.3$val),
                 sample = c(rep("before",nrow(interpolated.data.1)),rep("after_ln",nrow(interpolated.data.2)),
                            rep("after",nrow(interpolated.data.3))))
## scale data
df$value <- scale(df$value)
interpolated.data$val <- scale(interpolated.data$val)
summary(interpolated.data$val)
interpolated.data.1$val <- df[df$sample == "before","value"]
interpolated.data.2$val <- df[df$sample == "after_ln","value"]
interpolated.data.3$val <- df[df$sample == "after","value"]

## plot
p1 <- ggplot(interpolated.data.1, aes(x, 2e3 - y, color = val)) +
    geom_point(size = 0.1) +
    facet_wrap(~z, ncol = 1) +
    theme_void() +
    ggtitle("RPS9") + 
  #scale_color_viridis(discrete = F,limits = c(0,3))
  scale_color_gradientn(limits = c(-1,8),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))

p2 <- ggplot(interpolated.data.2, aes(x, 2e3 - y, color = val)) +
geom_point(size = 0.1) +
facet_wrap(~z, ncol = 1) +
theme_void() +
ggtitle("RPS9") + 
  #scale_color_viridis(discrete = F,limits = c(-1,8))
  #scale_color_gradientn(values = c(-1:1),colors = viridis(3,begin = 0, end = 0.5))
   scale_color_gradientn(limits = c(-1,8),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))

p3 <- ggplot(interpolated.data.3, aes(x, 2e3 - y, color = val)) +
geom_point(size = 0.1) +
facet_wrap(~z, ncol = 1) +
theme_void() +
ggtitle("RPS9") + 
  #scale_color_viridis(discrete = F,limits = c(0,3))
  #scale_color_gradientn(breaks = c(-Inf,1,3,5,Inf),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))
  scale_color_gradientn(limits = c(-1,8),colors = c("dark blue", "cyan", "yellow", "red", "dark red"))

pdf("results/RPS9_isoform_integrated.pdf")
print(p1)
print(p2)
print(p3)
dev.off()


## multiple 3d plots
p1 <- FeaturePlot3D(list_3d[[1]], features = "RPS9", scene = "scene", 
                    cols = c( "#3B0F70FF", "#8C2981FF","#DE4968FF"), 
                     pts.downsample = 5e4,
                    ,mode = "cloud", add.margins = 0,pt.size = 2)
p2 <- FeaturePlot3D(list_3d[[2]], features = "RPS9", scene = "scene2", 
                    cols = c("#3B0F70FF", "#8C2981FF","#DE4968FF","#FE9F6DFF"),
                    mode = "cloud", add.margins = 0,pt.size = 1.6)
p3 <- FeaturePlot3D(list_3d[[3]], features = "RPS9", scene = "scene2", 
                    cols = c("#3B0F70FF" ,"#DE4968FF" ,"#FE9F6DFF" ,"#FCFDBFFF"),
                    mode = "cloud", add.margins = 0,pt.size = 1.6)

plotly::subplot(p1,p2, p3,margin = 0)

