######### Libraries ##########
library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(biomaRt)
library(utils)
library(readr)
library(viridis)
library(ggridges)
library(ggalluvial)
library(ComplexHeatmap)
######## 1. read data #######
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_added_isoform.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_added_isoform.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_added_isoform.rds")

######## Spatial pattern ##########
### 1. Single isoform 2. Multi isoform (cover all clusters) 
### 3. Multi isoform (specific cluster) 4. Multi isoform (others)
###################################

####### 2 Spatial pattern classification  ######
spatial_pattern <- function(seurat,file_name){
  df <- as.data.frame(seurat@assays$ISO@counts)
  # add gene name
  df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
  t <- as.data.frame(table(df$geneId))
  # keep  1 to 0
  #df <- as.data.frame(seurat@assays$ISO@counts)
  df[df == 1] <- 0
  #  single isoform 
  single <- t[which(t$Freq==1),]
  df.single <- df[df$geneId %in% single$Var1,]
  df.single <- df.single[ , -which(names(df.single) %in% c("geneId"))]
  single_isoform <- rownames(df.single)
  dim(df.single)
  #  Multi isoform 
  
  # df <- df[rowSums(df) !=0,]
  # df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
  # 
  # delete num_zero >= 1% #spots
  # num_zero_df <- as.data.frame(rowSums(df == 0))
  # num_zero_df$isoform <- rownames(num_zero_df)
  # colnames(num_zero_df)[1] <- "num_zero"
  # num_zero_df <- num_zero_df[num_zero_df$num_zero <= 0.99*ncol(df) ,]
  # df <- df[num_zero_df$isoform,]
  # 
  multi <- t[which(t$Freq>1),]
  df.multi <- df[df$geneId %in% multi$Var1,]
  df.multi <- df.multi[ , -which(names(df.multi) %in% c("geneId"))]
  dim(df.multi)
  # 1031 492
  multi_all_isoform <- c()
  multi_one_isoform <- c()
  multi_others_isoform <- c()
  ## get cluster
  for (iso in rownames(df.multi)){
    cells <- t(df.multi[iso,])
    no_zero_cells <- names(cells[apply(cells, 1, function(row) all(row !=0 )), ])
    cell_table <- table(seurat@meta.data[no_zero_cells,"seurat_clusters"])
    ## filter <= 5
    #cell_table <- ifelse(cell_table <= 1,0,cell_table)
    if(!0 %in% cell_table){
      multi_all_isoform <- c(multi_all_isoform, iso)
    } else if (table(cell_table)["0"] == length(cell_table) - 1 ){
      multi_one_isoform <- c(multi_one_isoform, iso)
    } else{
      multi_others_isoform <- c(multi_others_isoform, iso)
    }
  }
  
  df_all_pattern <- data.frame(isoforms = c(single_isoform,multi_all_isoform,multi_one_isoform, multi_others_isoform), 
                               spatial_pattern = c(rep("single",length(single_isoform)),
                                                   rep("multi_all",length(multi_all_isoform)), 
                                                   rep("multi_one",length(multi_one_isoform)),
                                                   rep("multi_others",length(multi_others_isoform))))
  head(df_all_pattern)
  file_name <- paste0("results/",file_name,"_isoform_spatial_pattern_df.csv")
  write.table(df_all_pattern,file_name, sep=",",row.names = F)
  return(df_all_pattern)
}
isoform_before_sp <- spatial_pattern(before,"before_only1to0")
isoform_after_sp <- spatial_pattern(after,"after_only1to0")
isoform_after_ln_sp <- spatial_pattern(after_ln,"after_ln_only1to0")
####### 3 Spatial pattern stat and plot  （Fig 1B) ##########
### combine data 
isoform_before_sp <- as.data.frame(read_csv("results/isoform_total_stat/before_only1to0_isoform_spatial_pattern_df.csv"))
isoform_after_sp <- as.data.frame(read_csv("results/isoform_total_stat/after_only1to0_isoform_spatial_pattern_df.csv"))
isoform_after_ln_sp <- as.data.frame(read_csv("results/isoform_total_stat/after_ln_only1to0_isoform_spatial_pattern_df.csv"))

df_spatial_pattern_isoform <- rbind(isoform_before_sp,isoform_after_sp,isoform_after_ln_sp)
df_spatial_pattern_isoform$sample <- c(rep("before",nrow(isoform_before_sp)),rep("after",nrow(isoform_after_sp)),
                               rep("after_ln",nrow(isoform_after_ln_sp)))

dim(df_spatial_pattern_isoform)
table(df_spatial_pattern_isoform$spatial_pattern)
# multi_all    multi_one multi_others       single 
# 2445         1068        65023        18741 
write.table(df_spatial_pattern_isoform,"results/isoform_total_stat/spatial_pattern_isoform.csv", sep=",",row.names = F)
df_spatial_pattern_isoform <- as.data.frame(read_csv("results/isoform_total_stat/spatial_pattern_isoform.csv"))

### plot
pdf("results/spatial_pattern_distribution_isoform.pdf")
ggplot(df_spatial_pattern_isoform,aes(x = factor(sample,levels = c("before","after","after_ln")),
                              fill = spatial_pattern)) +  geom_histogram(stat = "count",position='fill') + 
  scale_fill_viridis(discrete = T) + xlab(label = "Frequency") + ylab("")
dev.off()

####### 4 # spots distribution  ###########
## function
num_spots_distribution <- function(seurat,seurat_name,isoforms_df){
  df <- as.data.frame(seurat@assays$ISO@counts)
  # delete <= 1
  #df[df == 1] <- 0
  #df <- df[rowSums(df) !=0,]
  
  isoforms_df <- isoforms_df[isoforms_df$sample == seurat_name,]
  df.multi.all <- df[rownames(df) %in% isoforms_df$isoforms,]
  #df.multi.all <- df.multi.all[ , -which(names(df.multi.all) %in% c("geneId"))]
  num_spots_df <- as.data.frame(Matrix::rowSums(df.multi.all))
  colnames(num_spots_df) <- c("num_spots")
  num_spots_df$cat <- ">1000"
  num_spots_df$cat[which(num_spots_df$num_spots<=1000)] <- "500-1000"
  num_spots_df$cat[which(num_spots_df$num_spots<500)] <- "200-500"
  num_spots_df$cat[which(num_spots_df$num_spots<200)] <- "100-200"
  num_spots_df$cat[which(num_spots_df$num_spots<100)] <- "50-100"
  num_spots_df$cat[which(num_spots_df$num_spots<50)] <- "25-50"
  num_spots_df$cat[which(num_spots_df$num_spots<25)] <- "10-25"
  num_spots_df$cat[which(num_spots_df$num_spots<10)] <- "<10"
  return(num_spots_df)
  
}
multi_one_num_spots_distribution <- function(seurat,seurat_name,isoforms_df){
  df <- as.data.frame(seurat@assays$ISO@counts)
  # delete <= 1
  #df[df == 1] <- 0
  #df <- df[rowSums(df) !=0,]
  isoforms_df <- isoforms_df[isoforms_df$sample == seurat_name,]
  df.multi.one <- df[rownames(df) %in% isoforms_df$isoforms,]
  #df.multi.one <- df.multi.one[ , -which(names(df.multi.one) %in% c("geneId"))]
  num_spots_df <- as.data.frame(Matrix::rowSums(df.multi.one))
  colnames(num_spots_df) <- c("num_spots")
  num_spots_df$cat <- ">500"
  num_spots_df$cat[which(num_spots_df$num_spots<=500)] <- "200-500"
  num_spots_df$cat[which(num_spots_df$num_spots<200)] <- "100-200"
  num_spots_df$cat[which(num_spots_df$num_spots<100)] <- "50-100"
  num_spots_df$cat[which(num_spots_df$num_spots<50)] <- "10-50"
  num_spots_df$cat[which(num_spots_df$num_spots<10)] <- "<10"
  return(num_spots_df)
  
}

## multi_all
multi_all_isoforms <- df_spatial_pattern_isoform[df_spatial_pattern_isoform$spatial_pattern == "multi_all",]
nrow(multi_all_isoforms)

before_num_spots_multi_all <- num_spots_distribution(before,"before",multi_all_isoforms)
after_num_spots_multi_all <- num_spots_distribution(after,"after",multi_all_isoforms)
after_ln_num_spots_multi_all <- num_spots_distribution(after_ln,"after_ln",multi_all_isoforms)
write.table(after_ln_num_spots_multi_all,"results/isoform_total_stat/after_ln_num_spots_multi_all.csv", sep=",",row.names = T)

## all
multi_all_all_num_spots <- rbind(before_num_spots_multi_all,after_num_spots_multi_all,after_ln_num_spots_multi_all)
multi_all_all_num_spots$sample <- c(rep("before",nrow(before_num_spots_multi_all)),
                                    rep("after",nrow(after_num_spots_multi_all)),
                              rep("after_ln",nrow(after_ln_num_spots_multi_all)))

## plot
pdf("results/num_spots_multi_all_all_samples.pdf", width = 15)
ggplot(multi_all_all_num_spots , aes(x = factor(cat,levels = c(">1000","500-1000","200-500","100-200",
    "50-100", "25-50", "10-25","<10")), fill =  factor(cat,levels = c(">1000","500-1000","200-500","100-200",
       "50-100", "25-50", "10-25","<10")))) + geom_bar() + xlab("Number of spots") +  theme_bw() + 
  theme(panel.grid = element_blank(),legend.position="none") +
  facet_wrap(factor(multi_all_all_num_spots$sample, levels = c("before","after","after_ln")), scales="free") +  scale_fill_brewer(palette = "RdPu",
                                                          direction = -1)
dev.off()

## multi_one
multi_one_isoforms <- df_spatial_pattern_isoform[df_spatial_pattern_isoform$spatial_pattern == "multi_one",]
nrow(multi_one_isoforms)

before_num_spots <- multi_one_num_spots_distribution(before,"before",multi_one_isoforms)
after_num_spots <- multi_one_num_spots_distribution(after,"after",multi_one_isoforms)
after_ln_num_spots <- multi_one_num_spots_distribution(after_ln,"after_ln",multi_one_isoforms)

## all
multi_one_all_num_spots <- rbind(before_num_spots,after_num_spots,after_ln_num_spots)
multi_one_all_num_spots$sample <- c(rep("before",nrow(before_num_spots)),rep("after",nrow(after_num_spots)),
                          rep("after_ln",nrow(after_ln_num_spots)))
## plot

pdf("results/num_spots_multi_one_all_samples.pdf", width = 15)
ggplot(multi_one_all_num_spots , aes(x = factor(cat,levels = c(">500","200-500","100-200",
        "50-100", "10-50","<10")), fill = factor(cat,levels = c(">500","200-500","100-200",
         "50-100", "10-50","<10")))) + geom_bar() + xlab("Number of spots") +  theme_bw() + 
           theme(panel.grid = element_blank(),legend.position="none") +
            facet_wrap(factor(multi_one_all_num_spots$sample, levels = c("before","after","after_ln")), scales="free") +  scale_fill_brewer(palette = "RdPu",
             direction = -1)
dev.off()


####### 5.  Average exp distribution  ###########
### multi_all
df.avg.exp <- df_spatial_pattern_isoform[df_spatial_pattern_isoform$spatial_pattern == "multi_all",]

avg_exp_distribution <- function(seurat,seurat_name,isoforms_df){
  df <- as.data.frame(seurat@assays$ISO@counts)
  #df[df == 1] <- 0
  isoforms_df <- isoforms_df[isoforms_df$sample == seurat_name,]
  df.avg.exp <- df[rownames(df) %in% isoforms_df$isoforms,]
  df.avg.exp[df.avg.exp == 0] <- NA
  df.avg.exp <- rowMeans(df.avg.exp,na.rm = T)
  return(df.avg.exp)
}

before_avg_exp_multi_all <- avg_exp_distribution(before,"before",multi_all_isoforms)
after_avg_exp_multi_all <- avg_exp_distribution(after,"after",multi_all_isoforms)
after_ln_avg_exp_multi_all <- avg_exp_distribution(after_ln,"after_ln",multi_all_isoforms)

## all
multi_all_all_avg_exp <- data.frame(avg_exp = c(before_avg_exp_multi_all,after_avg_exp_multi_all,after_ln_avg_exp_multi_all),
                                    sample = c(rep("before",length(before_avg_exp_multi_all)),
                                               rep("after",length(after_avg_exp_multi_all)),
                                          rep("after_ln",length(after_ln_avg_exp_multi_all))))
multi_all_all_avg_exp$avg_exp <- log10(multi_all_all_avg_exp$avg_exp)
summary(multi_all_all_avg_exp$avg_exp)
## plot
pdf("results/multi_all_avg_exp_all_samples_no1to0.pdf",height = 5)
ggplot(multi_all_all_avg_exp , aes(x = avg_exp, y = factor(sample, levels = c("before","after","after_ln")),fill = sample)) + 
  #geom_density() + scale_x_continuous(expand = c(0, 100)) + xlim(c(0,30))
  geom_density_ridges(scale = 2)   + xlab("log10(average experssion level)") + ylab("")
dev.off()

### multi_one
multi_one_isoforms <- df_spatial_pattern_isoform[df_spatial_pattern_isoform$spatial_pattern == "multi_one",]

before_avg_exp_multi_one <- avg_exp_distribution(before,"before",multi_one_isoforms)
after_avg_exp_multi_one <- avg_exp_distribution(after,"after",multi_one_isoforms)
after_ln_avg_exp_multi_one <- avg_exp_distribution(after_ln,"after_ln",multi_one_isoforms)

## all
multi_one_all_avg_exp <- data.frame(avg_exp = c(before_avg_exp_multi_one,after_avg_exp_multi_one,after_ln_avg_exp_multi_one),
                                    sample = c(rep("before",length(before_avg_exp_multi_one)),
                                               rep("after",length(after_avg_exp_multi_one)),
                                               rep("after_ln",length(after_ln_avg_exp_multi_one))))
multi_one_all_avg_exp$avg_exp <- log10(multi_one_all_avg_exp$avg_exp)
summary(multi_one_all_avg_exp$avg_exp)
## plot
pdf("results/multi_one_avg_exp_all_samples_no1to0.pdf",height = 5)
ggplot(multi_one_all_avg_exp , aes(x = avg_exp, y = factor(sample, levels = c("before","after","after_ln")),fill = sample)) + 
  #geom_density() + scale_x_continuous(expand = c(0, 100)) + xlim(c(0,30))
  geom_density_ridges(scale = 2)   + xlab("log10(average experssion level)") + ylab("")
dev.off()


####### 6. Check multi_all pattern ##########
multi_all_isoforms$nspots <- multi_all_all_num_spots$num_spots
multi_all_isoforms <- multi_all_isoforms[multi_all_isoforms$nspots >= 100,]
dim(multi_all_isoforms)
# 2445    3
####### 6.1 three samples overlapping and check avg exp level by heatmap #######
intersection_multi_all <- names(table(multi_all_isoforms$isoforms))[which(table(multi_all_isoforms$isoforms) == 3)]
intersection_multi_all <- multi_all_isoforms[multi_all_isoforms$isoforms %in% intersection_multi_all,]
dim(intersection_multi_all)
# 1356    3
table(intersection_multi_all$sample)
#after after_ln   before 
#452      452      452
## spots > 0.85
after_ln_num_spots_multi_all_num_spots_3000 <- after_ln_num_spots_multi_all[after_ln_num_spots_multi_all$num_spots >= quantile(after_ln_num_spots_multi_all$num_spots,0.85),]
after_ln_num_spots_multi_all_num_spots_3000 <- after_ln_num_spots_multi_all_num_spots_3000[rownames(after_ln_num_spots_multi_all_num_spots_3000) %in% intersection_multi_all$isoforms ,]
after_ln_num_spots_multi_all_num_spots_3000 <- after_ln_num_spots_multi_all_num_spots_3000[order(after_ln_num_spots_multi_all_num_spots_3000$num_spots,decreasing = T),]

after_num_spots_multi_all_num_spots_3000 <- after_num_spots_multi_all[after_num_spots_multi_all$num_spots >= quantile(after_num_spots_multi_all$num_spots,0.85),]
after_num_spots_multi_all_num_spots_3000 <- after_num_spots_multi_all_num_spots_3000[rownames(after_num_spots_multi_all_num_spots_3000) %in% intersection_multi_all$isoforms ,]
after_num_spots_multi_all_num_spots_3000 <- after_num_spots_multi_all_num_spots_3000[order(after_num_spots_multi_all_num_spots_3000$num_spots,decreasing = T),]

before_num_spots_multi_all_num_spots_3000 <- before_num_spots_multi_all[before_num_spots_multi_all$num_spots >= quantile(before_num_spots_multi_all$num_spots,0.85),]
before_num_spots_multi_all_num_spots_3000 <- before_num_spots_multi_all_num_spots_3000[rownames(before_num_spots_multi_all_num_spots_3000) %in% intersection_multi_all$isoforms ,]
before_num_spots_multi_all_num_spots_3000 <- before_num_spots_multi_all_num_spots_3000[order(before_num_spots_multi_all_num_spots_3000$num_spots,decreasing = T),]

intersection_multi_all_filtered <- intersect(rownames(after_ln_num_spots_multi_all_num_spots_3000),intersect(
                                             rownames(after_num_spots_multi_all_num_spots_3000),
                                             rownames(before_num_spots_multi_all_num_spots_3000)))

intersection_multi_all_filtered <- intersection_multi_all[intersection_multi_all$isoforms %in% intersection_multi_all_filtered, ]
dim(intersection_multi_all_filtered)
## function
get_exp_by_sample <- function(isoform_list,avg_exp_df,seurat,seurat_name){
  isoform_seurat <- isoform_list[isoform_list$sample == seurat_name,]
  df_seurat <- avg_exp_df[names(avg_exp_df) %in% isoform_seurat$isoforms]
  return(df_seurat)
}
## avg_level
before_count <- get_exp_by_sample(intersection_multi_all_filtered,before_avg_exp_multi_all,before,"before")
after_count <- get_exp_by_sample(intersection_multi_all_filtered,after_avg_exp_multi_all,after,"after")
after_ln_count <- get_exp_by_sample(intersection_multi_all_filtered,after_ln_avg_exp_multi_all,after_ln,"after_ln")

## combine
all_count <- data.frame(before = before_count, after = after_count, after_ln = after_ln_count)
typeof(all_count$before)
all_count$before <- as.numeric(all_count$before)
## scale
exp_all_count <- apply(all_count, 1, scale)
rownames(exp_all_count) <- colnames(all_count)
exp_all_count <- t(exp_all_count)
## kmeans 
km_exp_all_count = kmeans(exp_all_count, 4, nstart = 25)
km_exp_all_count <- km_exp_all_count$cluster
intersection_multi_all_filtered$cluster <- km_exp_all_count

metadata <- data.frame(exp = rownames(exp_all_count), cluster = km_exp_all_count )
metadata <- metadata[order(metadata$cluster,decreasing = F),]
metadata$cluster <- ifelse(metadata$cluster == 1,"Metastasis positive",
                           ifelse(metadata$cluster == 2,"Treatment targets",
                                  ifelse(metadata$cluster == 3,"Resistance","Metastasis negative")))
## heatmap (Fig 1C)
col <- brewer.pal(4, "PuBu")
names(col) <- c("Metastasis positive", "Treatment targets", "Resistance","Metastasis negative")
pdf("results/multi_all_overlapping_heatmap_0.85_no_row_names.pdf",height = 15,width = 10)
Heatmap(exp_all_count[,c("before","after","after_ln")], row_order = metadata$exp,
        , cluster_columns = F, row_km = 4, show_row_names = F,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
         labels = c("Metastasis positive", "Treatment targets", "Resistance","Metastasis negative"),
         labels_gp = gpar(col = "white", fontsize = 9)), col = col))
dev.off()

## spatial distribution
isoform_spatial_level <- function(seurat, gene_list, size ){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""), 
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=size, alpha = c(1, 1),ncol = 5)
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
isoform_spatial_level_plot <- function(gene){
  before_plot_list <- isoform_spatial_level(before, gene,2.5)
  after_plot_list <- isoform_spatial_level(after, gene,1.8)
  after_ln_plot_list <- isoform_spatial_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_all/overlapping/after_ln/isoform_spatial_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
    for(i in 1:length(before_plot_list)){
      print(before_plot_list[[i]])
    }}
  if (length(after_plot_list) != 0){
    for(i in 1:length(after_plot_list)){
      print(after_plot_list[[i]])
    }}
  if (length(after_ln_plot_list) != 0){
    for(i in 1:length(after_ln_plot_list)){
      print(after_ln_plot_list[[i]])
    }}
  dev.off()
}

multi_all_overlapping_after_ln <- c("RPL36","RPL31","RPL37A","RPS11","RPL7A","RPL24",
                                    "RPL27")
multi_all_overlapping_after <- c("RPS20","RPS3","RPL10","RPL12","RPS9")
for (gene in multi_all_overlapping_after){
  isoform_spatial_level_plot(gene)
}

####### 6.2 unique to some sample  #######
unique_multi_all <- names(table(multi_all_isoforms$isoforms))[which(table(multi_all_isoforms$isoforms) == 1)]
length(unique_multi_all)
unique_multi_all_df <- multi_all_isoforms[multi_all_isoforms$isoforms %in% unique_multi_all,]
table(unique_multi_all_df$sample)
# after after_ln   before 
# 38      466       99 

####### 6.3 unique to after ln and ##spots > 3000 ####
write.table(after_ln_unique_multi_all_df,"results/isoform_total_stat/after_ln_unique_multi_all_df.csv", sep=",",row.names = F)
after_ln_unique_multi_all_df <- unique_multi_all_df[unique_multi_all_df$sample == "after_ln",]
after_ln_num_spots_multi_all_num_spots_500 <- after_ln_num_spots_multi_all[after_ln_num_spots_multi_all$num_spots >= 3000,]
after_ln_num_spots_multi_all_num_spots_500 <- after_ln_num_spots_multi_all_num_spots_500[rownames(after_ln_num_spots_multi_all_num_spots_500) %in% after_ln_unique_multi_all_df$isoforms ,]
after_ln_num_spots_multi_all_num_spots_500 <- after_ln_num_spots_multi_all_num_spots_500[order(after_ln_num_spots_multi_all_num_spots_500$num_spots,decreasing = T),]
after_ln_num_spots_multi_all_num_spots_500$geneid <- sapply(strsplit(rownames(after_ln_num_spots_multi_all_num_spots_500), "\\.\\."), `[`, 1)
dim(after_ln_num_spots_multi_all_num_spots_500)
head(after_ln_num_spots_multi_all_num_spots_500)

isoform_spatial_level <- function(seurat, gene_list, size ){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""), 
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=size, alpha = c(1, 1),ncol = 5)
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
isoform_spatial_level_plot <- function(gene){
  before_plot_list <- isoform_spatial_level(before, gene,2.5)
  after_plot_list <- isoform_spatial_level(after, gene,1.8)
  after_ln_plot_list <- isoform_spatial_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_all/isoform_spatial_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
  for(i in 1:length(before_plot_list)){
    print(before_plot_list[[i]])
  }}
  if (length(after_plot_list) != 0){
  for(i in 1:length(after_plot_list)){
    print(after_plot_list[[i]])
  }}
  if (length(after_ln_plot_list) != 0){
  for(i in 1:length(after_ln_plot_list)){
    print(after_ln_plot_list[[i]])
  }}
  dev.off()
}

for (gene in after_ln_num_spots_multi_all_num_spots_500$geneid){
  isoform_spatial_level_plot(gene)
}
####### 6.4 unique to after  and ##spots > 0.8 ####
after_unique_multi_all_df <- unique_multi_all_df[unique_multi_all_df$sample == "after",]
after_num_spots_multi_all_num_spots_500 <- after_num_spots_multi_all[after_num_spots_multi_all$num_spots >= quantile(after_num_spots_multi_all$num_spots,0.8),]
after_num_spots_multi_all_num_spots_500 <- after_num_spots_multi_all_num_spots_500[rownames(after_num_spots_multi_all_num_spots_500) %in% after_unique_multi_all_df$isoforms ,]
after_num_spots_multi_all_num_spots_500 <- after_num_spots_multi_all_num_spots_500[order(after_num_spots_multi_all_num_spots_500$num_spots,decreasing = T),]
after_num_spots_multi_all_num_spots_500$geneid <- sapply(strsplit(rownames(after_num_spots_multi_all_num_spots_500), "\\.\\."), `[`, 1)
dim(after_num_spots_multi_all_num_spots_500)
head(after_num_spots_multi_all_num_spots_500)

isoform_spatial_level_plot <- function(gene){
  before_plot_list <- isoform_spatial_level(before, gene,2.5)
  after_plot_list <- isoform_spatial_level(after, gene,1.8)
  after_ln_plot_list <- isoform_spatial_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_all/unique_after/isoform_spatial_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
    for(i in 1:length(before_plot_list)){
      print(before_plot_list[[i]])
    }}
  if (length(after_plot_list) != 0){
    for(i in 1:length(after_plot_list)){
      print(after_plot_list[[i]])
    }}
  if (length(after_ln_plot_list) != 0){
    for(i in 1:length(after_ln_plot_list)){
      print(after_ln_plot_list[[i]])
    }}
  dev.off()
}
for (gene in after_num_spots_multi_all_num_spots_500$geneid){
  isoform_spatial_level_plot(gene)
}

####### 6.5 unique to before  and ##spots > 0.8 ####
before_unique_multi_all_df <- unique_multi_all_df[unique_multi_all_df$sample == "before",]
before_num_spots_multi_all_num_spots_500 <- before_num_spots_multi_all[before_num_spots_multi_all$num_spots >= quantile(before_num_spots_multi_all$num_spots,0.8),]
before_num_spots_multi_all_num_spots_500 <- before_num_spots_multi_all_num_spots_500[rownames(before_num_spots_multi_all_num_spots_500) %in% before_unique_multi_all_df$isoforms ,]
before_num_spots_multi_all_num_spots_500 <- before_num_spots_multi_all_num_spots_500[order(before_num_spots_multi_all_num_spots_500$num_spots,decreasing = T),]
before_num_spots_multi_all_num_spots_500$geneid <- sapply(strsplit(rownames(before_num_spots_multi_all_num_spots_500), "\\.\\."), `[`, 1)
dim(before_num_spots_multi_all_num_spots_500)
head(before_num_spots_multi_all_num_spots_500)

isoform_spatial_level_plot <- function(gene){
  before_plot_list <- isoform_spatial_level(before, gene,2.5)
  after_plot_list <- isoform_spatial_level(after, gene,1.8)
  after_ln_plot_list <- isoform_spatial_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_all/unique_before/isoform_spatial_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
    for(i in 1:length(before_plot_list)){
      print(before_plot_list[[i]])
    }}
  if (length(after_plot_list) != 0){
    for(i in 1:length(after_plot_list)){
      print(after_plot_list[[i]])
    }}
  if (length(after_ln_plot_list) != 0){
    for(i in 1:length(after_ln_plot_list)){
      print(after_ln_plot_list[[i]])
    }}
  dev.off()
}
for (gene in before_num_spots_multi_all_num_spots_500$geneid){
  isoform_spatial_level_plot(gene)
}



####### 6.6 RPS RPL family ##########
## add cluster info
intersection_multi_all_filtered$cluster <- metadata$cluster

# add family info
intersection_multi_all_filtered$geneid <- sapply(strsplit(intersection_multi_all_filtered$isoforms, "\\.\\."), `[`, 1)
intersection_multi_all_filtered$family <- ifelse(substr(x = intersection_multi_all_filtered$geneid, start = 1, stop = 3) == "RPL","RPL family",
                                             ifelse(substr(x = intersection_multi_all_filtered$geneid, start = 1, stop = 3) == "RPS","RPS family","Others"))
# Build link df
links <- intersection_multi_all_filtered[,c("cluster","family")]
links$num_spots <-  "Top 15% multi_all isoforms"
links <- links %>% group_by(cluster,family,num_spots) %>% summarize(Freq = n())
# colnames(links_1) <- c("source","target","value")
# links <- links_1
# 
# links$num_spots <-  rep("Top 15% multi_all isoforms",12)

# links <- rbind(links, data.frame(source = rep("Top 15% multi_all isoforms",4),
#                                             target = c("Metastasis negative",
#                                             "Metastasis positive", "Resistance",
#                                             "Treatment targets"), value = c(63,45,
#                                                                             66,21)))

#### ggalluvial (Fig 1D) #######
links_long <- to_lodes_form(data.frame(links),
                           key = "source",
                           axes = 1:3)

pdf("results/sankey_RPS_RPL_family.pdf",width = 9,height = 5)
ggplot(links, aes(axis1 = num_spots, axis2 = cluster, axis3 = family,  y = Freq)) + 
  geom_alluvium(aes(fill = cluster),width = 0.5,alpha = 0.5) +
  geom_stratum(aes(fill = after_stat(stratum)),width = 0.5) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
  theme_minimal() + 
  #scale_fill_brewer(type = "qual", palette = "Set2") + 
  scale_fill_viridis(discrete = T,begin = 0.2) +
  theme(legend.position = "none") + xlab("")+ylab ("Frequency")
dev.off()

####### 7. Check multi_one pattern ##########
multi_one_isoforms$nspots <- multi_one_all_num_spots$num_spots
multi_one_isoforms$geneid <- sapply(strsplit(multi_one_isoforms$isoforms, "\\.\\."), `[`, 1)
dim(multi_one_isoforms)
# 1068    4
## no overlapping 很正常
intersection_multi_one <- names(table(multi_one_isoforms$isoforms))[which(table(multi_one_isoforms$isoforms) == 3)]
intersection_multi_one <- multi_one_isoforms[multi_one_isoforms$isoforms %in% intersection_multi_one,]

dim(intersection_multi_one)
# 1356    3
## unique
table(multi_one_isoforms$sample)
#after after_ln   before 
#620      246      202 

## function
isoform_vlnplot_level <- function(seurat, gene_list, size ){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""), 
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = VlnPlot(seurat,group.by = "celltype", features = x,ncol = 5)
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
isoform_spatial_vlnplot <- function(gene){
  before_plot_list <- isoform_vlnplot_level(before, gene,2.5)
  after_plot_list <- isoform_vlnplot_level(after, gene,1.8)
  after_ln_plot_list <- isoform_vlnplot_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_one/isoform_vlnplot_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
    for(i in 1:length(before_plot_list)){
      print(before_plot_list[[i]])
    }}
  if (length(after_plot_list) != 0){
    for(i in 1:length(after_plot_list)){
      print(after_plot_list[[i]])
    }}
  if (length(after_ln_plot_list) != 0){
    for(i in 1:length(after_ln_plot_list)){
      print(after_ln_plot_list[[i]])
    }}
  dev.off()
}

isoform_spatial_level <- function(seurat, gene_list, size ){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""), 
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=size, alpha = c(1, 1),ncol = 5)
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
isoform_spatial_level_plot <- function(gene){
  before_plot_list <- isoform_spatial_level(before, gene,2.5)
  after_plot_list <- isoform_spatial_level(after, gene,1.8)
  after_ln_plot_list <- isoform_spatial_level(after_ln, gene,1.2)
  
  pdf(paste0("results/isoform_multi_one/isoform_spatial_three_samples_",gene,".pdf"), width=20, height=10)
  if (length(before_plot_list) != 0){
    for(i in 1:length(before_plot_list)){
      print(before_plot_list[[i]])
    }}
  if (length(after_plot_list) != 0){
    for(i in 1:length(after_plot_list)){
      print(after_plot_list[[i]])
    }}
  if (length(after_ln_plot_list) != 0){
    for(i in 1:length(after_ln_plot_list)){
      print(after_ln_plot_list[[i]])
    }}
  dev.off()
}
####### 7.0 cell annotation #######
### after LN 
after_ln@meta.data$celltype <- ifelse(after_ln$seurat_clusters == "0","Bcells, Myeloid and Tcells RGS1+",
                                      ifelse(after_ln$seurat_clusters == "1","Epithelial cells COL1A2+",
                                             ifelse(after_ln$seurat_clusters == "2","Bcells and Fibroblast IGLC2+",
                                                    ifelse(after_ln$seurat_clusters == "3","Bcells , Myeloid and Tcell APOC1+",
                                                           ifelse(after_ln$seurat_clusters == "4","Epithelial cells CXCL13+",
                                                                  ifelse(after_ln$seurat_clusters == "5","Epithelial cells FABP5+","oh no"))))))
### after 
after@meta.data$celltype <- ifelse(after$seurat_clusters == "0",'FRC and Epithelial COL1A2+',
                                   ifelse(after$seurat_clusters == "1",'Fibroblast, Pericytes and Endothelial MYH11+',
                                          ifelse(after$seurat_clusters == "2",'Bcell, Tcell, Myeloid IGHG4+',
                                                 ifelse(after$seurat_clusters == "3",'Fibroblast, Pericytes and Endothelial COL1A1-',
                                                        ifelse(after$seurat_clusters == "4",'Epithelial CRISP3+',
                                                               ifelse(after$seurat_clusters == "5",'FRC and Epithelial FABP5+',
                                                                      ifelse(after$seurat_clusters == "6",'Bcell, Tcell, Myeloid CXCL13+',
                                                                             ifelse(after$seurat_clusters == "7",'Epithelial MUC4+',"no"))))))))




table(after$seurat_clusters)
table(after$celltype)

## before
before@meta.data$celltype <- ifelse(before$seurat_clusters == "0",'FRC and Epithelial KRT178-',
                                    ifelse(before$seurat_clusters == "1",'FRC,Fibroblast and Epithelial CXCL14+',
                                           ifelse(before$seurat_clusters == "2",'FRC and Epithelial KRT4+',
                                                  ifelse(before$seurat_clusters == "3",'Bcell, FRC and Epithelial IFITM3+',
                                                         ifelse(before$seurat_clusters == "4",'FRC and Epithelial IGFBP5+',
                                                                ifelse(before$seurat_clusters == "5",'Epithelial MUC4+',"no"))))))


####### 7.1 unique to after_ln  #######
multi_one_after_ln <- multi_one_isoforms[multi_one_isoforms$sample == "after_ln",]
multi_one_after_ln <- multi_one_after_ln[order(multi_one_after_ln$nspots,decreasing = T),]
for (gene in multi_one_after_ln$geneid[1:3]){
  isoform_spatial_vlnplot(gene)
}
for (gene in multi_one_after_ln$geneid[1:10]){
  isoform_spatial_level_plot(gene)
}
multi_one_after_ln 
####### 7.2 unique to after  #######
multi_one_after <- multi_one_isoforms[multi_one_isoforms$sample == "after",]
multi_one_after <- multi_one_after[order(multi_one_after$nspots,decreasing = T),]
multi_one_after

####### 7.3 unique to before  #######
multi_one_before <- multi_one_isoforms[multi_one_isoforms$sample == "before",]
multi_one_before <- multi_one_before[order(multi_one_before$nspots,decreasing = T),]
multi_one_before



