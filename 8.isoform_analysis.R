######### Libraries ##########
library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(biomaRt)
###### 1. read seurat obj ######
## ori
#before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/3.Seurat/P1_before/P1_before.seurat.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after/p1_after.seurat.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after_LN/p1_after_LN.seurat.rds")

## added info
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_added_isoform.rds")
##### 2. cell annotation #######
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


###### 3. ISOG assay containing the gene-level nanopore data  #####
## before
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/before_minion/Analysis/ont_A_10_28_genematrix.txt", stringsAsFactors = F)
data = all[,2:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) <- all$geneId
vect <- colnames(before@assays$Spatial)
data <- data[vect]
before[["ISOG"]] <- CreateAssayObject(counts = data)
before <- NormalizeData(object = before, assay = "ISOG")
before <- ScaleData(object = before, assay = "ISOG")

## after
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/after_minion/after_genematrix.txt", stringsAsFactors = F)
data = all[,2:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) <- all$geneId
vect <- colnames(after@assays$Spatial)
data <- data[vect]
after[["ISOG"]] <- CreateAssayObject(counts = data)
after <- NormalizeData(object = after, assay = "ISOG")
after <- ScaleData(object = after, assay = "ISOG")

## afte-LN
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/after_minion/after_LN_genematrix.txt", stringsAsFactors = F)
data = all[,2:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) <- all$geneId
vect <- colnames(after_ln@assays$Spatial)
data <- data[vect]
after_ln[["ISOG"]] <- CreateAssayObject(counts = data)
after_ln <- NormalizeData(object = after_ln, assay = "ISOG")
after_ln <- ScaleData(object = after_ln, assay = "ISOG")

####### 4. ISO assay containing the isoform-level transcript information (only molecules where all eons are ovserved are kept) ########
## before
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/before_minion/Analysis/ont_A_10_28_isomatrix.txt", stringsAsFactors = F)
data = all[,3:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) = paste(all$geneId, all$transcriptId, sep="..")
vect <- colnames(before@assays$Spatial)
idx <- grep("undef", rownames(data), invert = TRUE)
data <- data[idx,vect]
dim(data)
before[["ISO"]] <- CreateAssayObject(counts = data)
before <- NormalizeData(object = before, assay = "ISO")
before <- ScaleData(before, assay = "ISO")
## save seurat object afted adding isoform information 
saveRDS(before,"~/onedrive/Work/phD/phd_project/SiT/results/before_added_isoform.rds")

## after
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/after_minion/after_isomatrix.txt", stringsAsFactors = F)
data = all[,3:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) = paste(all$geneId, all$transcriptId, sep="..")
vect <- colnames(after@assays$Spatial)
idx <- grep("undef", rownames(data), invert = TRUE)
data <- data[idx,vect]
dim(data)
after[["ISO"]] <- CreateAssayObject(counts = data)
after <- NormalizeData(object = after, assay = "ISO")
after <- ScaleData(after, assay = "ISO")
## save seurat object afted adding isoform information 
saveRDS(after,"~/onedrive/Work/phD/phd_project/SiT/results/after_added_isoform.rds")

## after_ln
all = read.delim("~/onedrive/Work/phD/phd_project/SiT/rawData/after_minion/after_ln_isomatrix.txt", stringsAsFactors = F)
data = all[,3:ncol(all)]
colnames(data) <- paste(colnames(data),"-1", sep="")
rownames(data) = paste(all$geneId, all$transcriptId, sep="..")
vect <- colnames(after_ln@assays$Spatial)
idx <- grep("undef", rownames(data), invert = TRUE)
data <- data[idx,vect]
dim(data)
after_ln[["ISO"]] <- CreateAssayObject(counts = data)
after_ln <- NormalizeData(object = after_ln, assay = "ISO")
after_ln <- ScaleData(after_ln, assay = "ISO")
## save seurat object afted adding isoform information 
saveRDS(after_ln,"~/onedrive/Work/phD/phd_project/SiT/results/after_ln_added_isoform.rds")


################### Isoform Spatially distribution stat  ###################
df <- as.data.frame(before@assays$ISO@counts)
df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
df$transcriptId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 2)
df <- as.data.frame(table(sapply(strsplit(rownames(before@assays$ISO@counts), "\\.\\."), `[`, 1)))
t <- as.data.frame(table(df$Freq))
t$rep <- "before"

t$pct <- round(100*t$Freq/sum(t$Freq), digits = 1)
dat <- t
pdf("results/before_minion/isoform.per.gene.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(data=dat, aes(x=Var1, y=Freq, fill=rep)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Freq), vjust=-0.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()
dev.off()

######## 5. generate MULTI assay containing only isoforms from multi-isoforms genes######
#### before ####
df <- as.data.frame(before@assays$ISO@counts)
df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df$geneId))
multi <- t[which(t$Freq>1),]
df.multi <- df[df$geneId %in% multi$Var1,]
df.multi <- df.multi[ , -which(names(df.multi) %in% c("geneId"))]
somme <- as.data.frame(Matrix::rowSums(df.multi))
colnames(somme) <- c("somme")
somme$cat <- ">1000"
somme$cat[which(somme$somme<=1000)] <- "500-1000"
somme$cat[which(somme$somme<500)] <- "200-500"
somme$cat[which(somme$somme<200)] <- "100-200"
somme$cat[which(somme$somme<100)] <- "50-100"
somme$cat[which(somme$somme<50)] <- "25-50"
somme$cat[which(somme$somme<25)] <- "10-25"
somme$cat[which(somme$somme<10)] <- "<10"
#new_order <- with(dat, reorder(variety , note, median , na.rm=T))
pdf("results/before_minion/multi_isoforms_number.per.gene.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(somme, aes(cat)) +
  geom_bar(fill = "steelblue")
dev.off()
##
one_thousand_somme <- somme[somme$cat == "500-1000",]
one_thousand_somme <- one_thousand_somme[order(one_thousand_somme$somme,decreasing = T),]
one_thousand_somme <- rownames(one_thousand_somme)
### Majority isoform 
df.multi.high <- df.multi
df.multi.high$geneId <- sapply(strsplit(rownames(df.multi.high), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df.multi.high$geneId))
multi <- t[which(t$Freq>1),]
df.multi.high <- df.multi.high[df.multi.high$geneId %in% multi$Var1,]
df.multi.high <- df.multi.high[ , -which(names(df.multi.high) %in% c("geneId"))]
gg <- multi$Var1
majoritary <- character()
dat <- data.frame()
for(i in 1:length(gg)){
  idx <- grep(paste0("^",gg[i],"\\.\\."),rownames(df.multi.high))
  iso <- Matrix::rowSums(df.multi.high[idx,])
  x <- max(iso)/sum(iso)
  
  majoritary <- c(majoritary, names(which.max(iso)))
  
  de <- data.frame(gg[i],max(iso),sum(iso),x)
  names(de)<-c("gene","max","total","ratiomax")
  dat <- rbind(dat, de)
  
  #print(paste0(gg[i], ",", max(iso), ",", sum(iso), ",", x, ",", sep=""))
}
hs=hist(dat$ratiomax, breaks=20)
pdf("results/before_minion/majority_isoform_ratio.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(data=dat, aes(x=total,y=ratiomax)) +  
  geom_point(shape = 21, colour = "black", fill="lightsteelblue3") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  scale_x_continuous(trans = 'log10') +
  geom_text_repel(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), aes(label=gene), size=5) +
  geom_point(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), col="red") 
dev.off()
## add info
before[["MULTI"]] <- CreateAssayObject(counts = df.multi.high)
before <- NormalizeData(object = before, assay = "MULTI")
before <- ScaleData(before, assay="MULTI")
##save

### Isoform switch detection
DefaultAssay(object = before) <- "MULTI"
total <- data.frame()
totaladj <- data.frame()
clusters <- unique(before@meta.data$celltype)

for (i in (1:5)){
  k <- i+1
  for (j in (k:6)){
    if(i != j){
      print(paste(i, " ", j, " ",clusters[i], " vs ", clusters[j], sep=""))
      
      markers <- FindMarkers(object = before, ident.1=clusters[i], ident.2=clusters[j],group.by =  "celltype")
      markers$cluster <- clusters[j]
      markers$region.1 <- clusters[i]
      markers$region.2 <- clusters[j]
      #markers$contrast <- paste(clusters[i], "vs", clusters[j], sep=" ")
      markers[which(markers$avg_log2FC>0),]$cluster <- clusters[i]
      markers$geneId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 1)
      markers$transcriptId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 2)
      markers$majo <- 0
      markers$majo[which(rownames(markers) %in% majoritary)] <- 1
      markers <- subset(markers, markers$majo == 1 | (markers$p_val < 0.05 & markers$majo == 0))
      all.genes <- unique(markers$geneId)
      for (k in (1:length(all.genes))){
        sub <- markers[which(markers$geneId == all.genes[k]),]
        nb.clusters <- unique(sub$cluster)
        nb.transcripts <- unique(sub$transcriptId)
        
        if(length(nb.clusters) > 1 & length(nb.transcripts) > 1){
          totaladj <- rbind(totaladj, sub)
        }
      }
      print (dim(totaladj))
      print (unique(totaladj$geneId))
    }
  }
}
write.table(totaladj, file="results/before_minion/before.isoswitch.csv", sep=",",row.names = F)

totaladj$abs_log2FC <- abs(totaladj$avg_log2FC)
totaladj <- totaladj[order(totaladj$abs_log2FC,decreasing = T),]

### really sig
totaladj_sig <- totaladj[totaladj$p_val_adj <= 0.05,]
totaladj_sig <- unique(totaladj_sig$geneId)
#### after ####
df <- as.data.frame(after@assays$ISO@counts)
df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df$geneId))
multi <- t[which(t$Freq>1),]
df.multi <- df[df$geneId %in% multi$Var1,]
df.multi <- df.multi[ , -which(names(df.multi) %in% c("geneId"))]
somme <- as.data.frame(Matrix::rowSums(df.multi))
colnames(somme) <- c("somme")
somme$cat <- ">1000"
somme$cat[which(somme$somme<=1000)] <- "500-1000"
somme$cat[which(somme$somme<500)] <- "200-500"
somme$cat[which(somme$somme<200)] <- "100-200"
somme$cat[which(somme$somme<100)] <- "50-100"
somme$cat[which(somme$somme<50)] <- "25-50"
somme$cat[which(somme$somme<25)] <- "10-25"
somme$cat[which(somme$somme<10)] <- "<10"
#new_order <- with(dat, reorder(variety , note, median , na.rm=T))
pdf("results/after_minion/multi_isoforms_number.per.gene.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(somme, aes(cat)) +
  geom_bar(fill = "steelblue")
dev.off()
##
one_thousand_somme <- somme[somme$cat == "500-1000",]
one_thousand_somme <- one_thousand_somme[order(one_thousand_somme$somme,decreasing = T),]
one_thousand_somme <- rownames(one_thousand_somme)
### Majority isoform 
df.multi.high <- df.multi
df.multi.high$geneId <- sapply(strsplit(rownames(df.multi.high), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df.multi.high$geneId))
multi <- t[which(t$Freq>1),]
df.multi.high <- df.multi.high[df.multi.high$geneId %in% multi$Var1,]
df.multi.high <- df.multi.high[ , -which(names(df.multi.high) %in% c("geneId"))]
gg <- multi$Var1
majoritary <- character()
dat <- data.frame()
for(i in 1:length(gg)){
  idx <- grep(paste0("^",gg[i],"\\.\\."),rownames(df.multi.high))
  iso <- Matrix::rowSums(df.multi.high[idx,])
  x <- max(iso)/sum(iso)
  
  majoritary <- c(majoritary, names(which.max(iso)))
  
  de <- data.frame(gg[i],max(iso),sum(iso),x)
  names(de)<-c("gene","max","total","ratiomax")
  dat <- rbind(dat, de)
  
  #print(paste0(gg[i], ",", max(iso), ",", sum(iso), ",", x, ",", sep=""))
}
hs=hist(dat$ratiomax, breaks=20)
pdf("results/after_minion/majority_isoform_ratio.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(data=dat, aes(x=total,y=ratiomax)) +  
  geom_point(shape = 21, colour = "black", fill="lightsteelblue3") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  scale_x_continuous(trans = 'log10') +
  geom_text_repel(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), aes(label=gene), size=5) +
  geom_point(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), col="red") 
dev.off()
## add info
after[["MULTI"]] <- CreateAssayObject(counts = df.multi.high)
after <- NormalizeData(object = after, assay = "MULTI")
after <- ScaleData(after, assay="MULTI")
##save
saveRDS(after,"~/onedrive/Work/phD/phd_project/SiT/results/after_added_isoform.rds")

### Isoform switch detection
DefaultAssay(object = after) <- "MULTI"
total <- data.frame()
totaladj <- data.frame()
clusters <- unique(after@meta.data$celltype)

for (i in (1:5)){
  k <- i+1
  for (j in (k:6)){
    if(i != j){
      print(paste(i, " ", j, " ",clusters[i], " vs ", clusters[j], sep=""))
      
      markers <- FindMarkers(object = after, ident.1=clusters[i], ident.2=clusters[j],group.by =  "celltype")
      markers$cluster <- clusters[j]
      markers$region.1 <- clusters[i]
      markers$region.2 <- clusters[j]
      #markers$contrast <- paste(clusters[i], "vs", clusters[j], sep=" ")
      markers[which(markers$avg_log2FC>0),]$cluster <- clusters[i]
      markers$geneId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 1)
      markers$transcriptId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 2)
      markers$majo <- 0
      markers$majo[which(rownames(markers) %in% majoritary)] <- 1
      markers <- subset(markers, markers$majo == 1 | (markers$p_val < 0.05 & markers$majo == 0))
      all.genes <- unique(markers$geneId)
      for (k in (1:length(all.genes))){
        sub <- markers[which(markers$geneId == all.genes[k]),]
        nb.clusters <- unique(sub$cluster)
        nb.transcripts <- unique(sub$transcriptId)
        
        if(length(nb.clusters) > 1 & length(nb.transcripts) > 1){
          totaladj <- rbind(totaladj, sub)
        }
      }
      print (dim(totaladj))
      print (unique(totaladj$geneId))
    }
  }
}
write.table(totaladj, file="results/after_minion/after.isoswitch.csv", sep=",",row.names = F)

# totaladj$abs_log2FC <- abs(totaladj$avg_log2FC)
# totaladj <- totaladj[order(totaladj$abs_log2FC,decreasing = T),]
# 
# ### really sig
# totaladj_sig <- totaladj[totaladj$p_val_adj <= 0.05,]
# totaladj_sig <- unique(totaladj_sig$geneId)
#### afterLN ####
df <- as.data.frame(after_ln@assays$ISO@counts)
df$geneId <- sapply(strsplit(rownames(df), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df$geneId))
multi <- t[which(t$Freq>1),]
df.multi <- df[df$geneId %in% multi$Var1,]
df.multi <- df.multi[ , -which(names(df.multi) %in% c("geneId"))]
somme <- as.data.frame(Matrix::rowSums(df.multi))
colnames(somme) <- c("somme")
somme$cat <- ">1000"
somme$cat[which(somme$somme<=1000)] <- "500-1000"
somme$cat[which(somme$somme<500)] <- "200-500"
somme$cat[which(somme$somme<200)] <- "100-200"
somme$cat[which(somme$somme<100)] <- "50-100"
somme$cat[which(somme$somme<50)] <- "25-50"
somme$cat[which(somme$somme<25)] <- "10-25"
somme$cat[which(somme$somme<10)] <- "<10"
#new_order <- with(dat, reorder(variety , note, median , na.rm=T))
pdf("results/after_minion/after_ln_multi_isoforms_number.per.gene.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(somme, aes(cat)) +
  geom_bar(fill = "steelblue")
dev.off()
##
one_thousand_somme <- somme[somme$cat == "500-1000",]
one_thousand_somme <- one_thousand_somme[order(one_thousand_somme$somme,decreasing = T),]
one_thousand_somme <- rownames(one_thousand_somme)
### Majority isoform 
df.multi.high <- df.multi
df.multi.high$geneId <- sapply(strsplit(rownames(df.multi.high), "\\.\\."), `[`, 1)
t <- as.data.frame(table(df.multi.high$geneId))
multi <- t[which(t$Freq>1),]
df.multi.high <- df.multi.high[df.multi.high$geneId %in% multi$Var1,]
df.multi.high <- df.multi.high[ , -which(names(df.multi.high) %in% c("geneId"))]
gg <- multi$Var1
majoritary <- character()
dat <- data.frame()
for(i in 1:length(gg)){
  idx <- grep(paste0("^",gg[i],"\\.\\."),rownames(df.multi.high))
  iso <- Matrix::rowSums(df.multi.high[idx,])
  x <- max(iso)/sum(iso)
  
  majoritary <- c(majoritary, names(which.max(iso)))
  
  de <- data.frame(gg[i],max(iso),sum(iso),x)
  names(de)<-c("gene","max","total","ratiomax")
  dat <- rbind(dat, de)
  
  #print(paste0(gg[i], ",", max(iso), ",", sum(iso), ",", x, ",", sep=""))
}
hs=hist(dat$ratiomax, breaks=20)
pdf("results/after_minion/after_ln_majority_isoform_ratio.pdf", width=8, height=6, useDingbats=FALSE)
ggplot(data=dat, aes(x=total,y=ratiomax)) +  
  geom_point(shape = 21, colour = "black", fill="lightsteelblue3") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  scale_x_continuous(trans = 'log10') +
  geom_text_repel(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), aes(label=gene), size=5) +
  geom_point(data=subset(dat, dat$ratiomax < 0.5 & dat$total > 1000), col="red") 
dev.off()
## add info
after_ln[["MULTI"]] <- CreateAssayObject(counts = df.multi.high)
after_ln <- NormalizeData(object = after_ln, assay = "MULTI")
after_ln <- ScaleData(after_ln, assay="MULTI")
##save
saveRDS(after_ln,"~/onedrive/Work/phD/phd_project/SiT/results/after_ln_added_isoform.rds")

### Isoform switch detection
DefaultAssay(object = after_ln) <- "MULTI"
total <- data.frame()
totaladj <- data.frame()
clusters <- unique(after_ln@meta.data$celltype)

for (i in (1:5)){
  k <- i+1
  for (j in (k:6)){
    if(i != j){
      print(paste(i, " ", j, " ",clusters[i], " vs ", clusters[j], sep=""))
      
      markers <- FindMarkers(object = after_ln, ident.1=clusters[i], ident.2=clusters[j],group.by =  "celltype")
      markers$cluster <- clusters[j]
      markers$region.1 <- clusters[i]
      markers$region.2 <- clusters[j]
      #markers$contrast <- paste(clusters[i], "vs", clusters[j], sep=" ")
      markers[which(markers$avg_log2FC>0),]$cluster <- clusters[i]
      markers$geneId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 1)
      markers$transcriptId <- sapply(strsplit(rownames(markers), "\\.\\."), `[`, 2)
      markers$majo <- 0
      markers$majo[which(rownames(markers) %in% majoritary)] <- 1
      markers <- subset(markers, markers$majo == 1 | (markers$p_val < 0.05 & markers$majo == 0))
      all.genes <- unique(markers$geneId)
      for (k in (1:length(all.genes))){
        sub <- markers[which(markers$geneId == all.genes[k]),]
        nb.clusters <- unique(sub$cluster)
        nb.transcripts <- unique(sub$transcriptId)
        
        if(length(nb.clusters) > 1 & length(nb.transcripts) > 1){
          totaladj <- rbind(totaladj, sub)
        }
      }
      print (dim(totaladj))
      print (unique(totaladj$geneId))
    }
  }
}
write.table(totaladj, file="results/after_minion/after_ln.isoswitch.csv", sep=",",row.names = F)
totaladj <- read_csv("results/after_minion/after_ln.isoswitch.csv")
# totaladj$abs_log2FC <- abs(totaladj$avg_log2FC)
# totaladj <- totaladj[order(totaladj$abs_log2FC,decreasing = T),]
# 
# ### really sig
# totaladj_sig <- totaladj[totaladj$p_val_adj <= 0.05,]
# totaladj_sig <- unique(totaladj_sig$geneId)


########## 6. Check isoform number vs cluster ######
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/before_minion/before_added_isoform.rds")
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_minion/after_added_isoform.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_added_isoform.rds")

##### before 
df <- as.data.frame(before@assays$MULTI@counts)
before_isoform_number <- c()
for (cluster in unique(before$seurat_clusters)){
  cluster_isoform <- df[,rownames(before@meta.data[before$seurat_clusters == cluster,])]
  cluster_isoform <- cluster_isoform[rowSums(cluster_isoform) > 0,]
  before_isoform_number <- c(before_isoform_number,nrow(cluster_isoform))
}
before_isoform_number
#ISO 10190 10191  6446  7090 12567  7708
#MULTI 8038 8067 5164 5528 9818 6149
##### after 
df <- as.data.frame(after@assays$MULTI@counts)
after_isoform_number <- c()
for (cluster in unique(after$seurat_clusters)){
  cluster_isoform <- df[,rownames(after@meta.data[after$seurat_clusters == cluster,])]
  cluster_isoform <- cluster_isoform[rowSums(cluster_isoform) > 0,]
  after_isoform_number <- c(after_isoform_number,nrow(cluster_isoform))
}
after_isoform_number
#ISO 18768  8677 12415  7608 10835 22950 18458 12983
#MULTi 15662  7226 10255  6412  9136 18889 15089 10806
##### after_ln 
df <- as.data.frame(after_ln@assays$MULTI@counts)
after_ln_isoform_number <- c()
for (cluster in unique(after_ln$seurat_clusters)){
  cluster_isoform <- df[,rownames(after_ln@meta.data[after_ln$seurat_clusters == cluster,])]
  cluster_isoform <- cluster_isoform[rowSums(cluster_isoform) > 0,]
  after_ln_isoform_number <- c(after_ln_isoform_number,nrow(cluster_isoform))
}
after_ln_isoform_number
#ISO 17325 22655 16038 21141 15899 16560
#MULTI 14179 18485 13238 17256 13128 13673

########## 7. Check isoform marker for each cluster #######
##### after_ln
after_ln_added_isoform <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/after_ln_added_isoform.rds")

cluster_markers <- FindAllMarkers(after_ln_added_isoform,assay = "MULTI")
cluster_markers_sig <- cluster_markers[cluster_markers$p_val <= 0.05,]
dim(cluster_markers_sig)
top50 <- cluster_markers_sig %>% group_by(cluster) %>% top_n(50, avg_log2FC)

write.csv(cluster_markers,"~/onedrive/Work/phD/phd_project/SiT/results/isoform_cluster_markers/after_ln_cluster_markers.csv")
write.csv(cluster_markers_sig,"~/onedrive/Work/phD/phd_project/SiT/results/isoform_cluster_markers/after_ln_cluster_markers_sig.csv")




########## 8. check isoform spatial level #######
## function
isoform_spatial_level <- function(gene_list,seurat){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""), rownames(seurat@assays$ISO@data))])
    DefaultAssay(before) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.3, alpha = c(1, 1),ncol = 3)
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
########### CD74 only in after ISO marker #####
plot_list <- list()
seurat <- after
for (i in 1:length(gene_list)){
  x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^","CD74","..",sep=""),
                                                    rownames(seurat@assays$ISO@data))])
  DefaultAssay(seurat) <- "ISO"
  if (length (x) != 0) {
    p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.8, alpha = c(1, 1),ncol = )
    plot_list[[i]] <- p
  }}

filename <- "~/onedrive/Work/phD/phd_project/SiT/results/after_minion/isoform_marker_after_ln_CD74_spatial_distribution.pdf"
pdf(filename, width=20, height=10)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()

########### CD74 only in before ISO marker #####
plot_list <- list()
seurat <- before
for (i in 1:length(gene_list)){
  x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^","CD74","..",sep=""),
                                                    rownames(seurat@assays$ISO@data))])
  DefaultAssay(seurat) <- "ISO"
  if (length (x) != 0) {
    p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=2.8, alpha = c(1, 1),ncol = )
    plot_list[[i]] <- p
  }}

filename <- "~/onedrive/Work/phD/phd_project/SiT/results/before_minion/isoform_marker_after_ln_CD74_spatial_distribution.pdf"
pdf(filename, width=20, height=10)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()


######## RPS/L isoform overlap in CD8_CXCL13 clusters ########
### after_ln RPL
after_ln_cluster_markers <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/isoform_cluster_markers/after_ln_cluster_markers.csv")
after_ln_cluster_markers <- after_ln_cluster_markers[after_ln_cluster_markers$cluster == 0,]
after_ln_cluster_markers <- after_ln_cluster_markers[after_ln_cluster_markers$p_val_adj <= 0.05,]
after_ln_cluster_markers_RPL <- after_ln_cluster_markers[grep("RPL..", after_ln_cluster_markers$X),]
summary(abs(after_ln_cluster_markers_RPL$avg_log2FC))

gene_list <- unique(sapply(strsplit(after_ln_cluster_markers_RPL$gene, "\\.\\."), `[`, 1))

plot_list <- list()
seurat <- after_ln
for (i in 1:length(gene_list)){
  x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0(gene_list[i],"\\.\\.",sep=""),
                                                    rownames(seurat@assays$ISO@data))])
  DefaultAssay(seurat) <- "ISO"
  if (length (x) != 0) {
    p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.3, alpha = c(1, 1),ncol = )
    plot_list[[i]] <- p
  }}


filename <- "~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/cd8_cxcl13_isoform_RPL_spatial_distribution.pdf"
pdf(filename, width=20, height=10)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()
# RPL5 !!
isoform_spatial_level <- function(gene_list,seurat){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""),
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=2.3, alpha = c(1, 1),ncol = )
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
isoform_vlnplot_level <- function(gene_list,seurat){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""),
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = VlnPlot(seurat, features = x ,group.by = "celltype")
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
filename <- "~/onedrive/Work/phD/phd_project/SiT/results/isoform_test/RPL5_spatial_distribution.pdf"
plot_list_after_ln <- isoform_spatial_level(c("RPL5"),after_ln)
plot_list_after <- isoform_spatial_level(c("RPL5"),after)
plot_list_before <- isoform_spatial_level(c("RPL5"),before)

pdf(filename, width=20, height=10)
for(i in 1:length(plot_list_before)){
  print(plot_list_before[[i]])
}
for(i in 1:length(plot_list_after)){
  print(plot_list_after[[i]])
  
}
for(i in 1:length(plot_list_after_ln)){
  print(plot_list_after_ln[[i]])
}

dev.off()

## dotplot and Vlnplot
filename <- "~/onedrive/Work/phD/phd_project/SiT/results/isoform_test/RPL5_vlnplot_distribution.pdf"
plot_list_after_ln <- isoform_vlnplot_level(c("RPL5"),after_ln)
plot_list_after <- isoform_vlnplot_level(c("RPL5"),after)
plot_list_before <- isoform_vlnplot_level(c("RPL5"),before)

pdf(filename, width=20, height=50)
for(i in 1:length(plot_list_before)){
  print(plot_list_before[[i]])
}
for(i in 1:length(plot_list_after)){
  print(plot_list_after[[i]])
  
}
for(i in 1:length(plot_list_after_ln)){
  print(plot_list_after_ln[[i]])
}

dev.off()


# RPL11 !!
isoform_spatial_level <- function(gene_list,seurat){
  plot_list <- list()
  for (i in 1:length(gene_list)){
    x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0("^",gene_list[i],"..",sep=""),
                                                      rownames(seurat@assays$ISO@data))])
    DefaultAssay(seurat) <- "ISO"
    if (length (x) != 0) {
      p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.3, alpha = c(1, 1),ncol = )
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}
filename <- "~/onedrive/Work/phD/phd_project/SiT/results/isoform_test/RPL11_spatial_distribution.pdf"
plot_list_after_ln <- isoform_spatial_level(c("RPL11"),after_ln)
plot_list_after <- isoform_spatial_level(c("RPL11"),after)
plot_list_before <- isoform_spatial_level(c("RPL11"),before)

pdf(filename, width=20, height=10)
for(i in 1:length(plot_list_before)){
  print(plot_list_before[[i]])
}
for(i in 1:length(plot_list_after)){
  print(plot_list_after[[i]])
  
}
for(i in 1:length(plot_list_after_ln)){
  print(plot_list_after_ln[[i]])
}

dev.off()


pdf("~/onedrive/Work/phD/phd_project/SiT/results/RPL5_gene_spatial_distribution.pdf")
SpatialFeaturePlot(before, features = "RPL5", pt.size.factor=2.3, alpha = c(1, 1))
SpatialFeaturePlot(after, features = "RPL5", pt.size.factor=1.8, alpha = c(1, 1))
SpatialFeaturePlot(after_ln, features = "RPL5", pt.size.factor=1.3, alpha = c(1, 1))
dev.off()


######### Malignant epi cell marker isoform & RPL/S

after_ln_cluster_markers <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/isoform_cluster_markers/after_ln_cluster_markers.csv")
after_ln_cluster_markers <- after_ln_cluster_markers[after_ln_cluster_markers$cluster == 5,]
after_ln_cluster_markers <- after_ln_cluster_markers[after_ln_cluster_markers$p_val_adj <= 0.05,]
after_ln_cluster_markers_RPL <- after_ln_cluster_markers[grep("RPL..", after_ln_cluster_markers$X),]
summary(abs(after_ln_cluster_markers_RPL$avg_log2FC))

gene_list <- unique(sapply(strsplit(after_ln_cluster_markers_RPL$gene, "\\.\\."), `[`, 1))

plot_list <- list()
seurat <- after_ln
for (i in 1:length(gene_list)){
  x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0(gene_list[i],"\\.\\.",sep=""),
                                                    rownames(seurat@assays$ISO@data))])
  DefaultAssay(seurat) <- "ISO"
  if (length (x) != 0) {
    p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.3, alpha = c(1, 1),ncol = )
    plot_list[[i]] <- p
  }}


filename <- "~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/malignant_isoform_RPL_spatial_distribution.pdf"
pdf(filename, width=20, height=10)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()

after_ln_cluster_markers_RPS <- after_ln_cluster_markers[grep("RPS..", after_ln_cluster_markers$X),]
summary(abs(after_ln_cluster_markers_RPS$avg_log2FC))

gene_list <- unique(sapply(strsplit(after_ln_cluster_markers_RPS$gene, "\\.\\."), `[`, 1))

plot_list <- list()
seurat <- after_ln
for (i in 1:length(gene_list)){
  x <- unique(rownames(seurat@assays$ISO@data)[grep(paste0(gene_list[i],"\\.\\.",sep=""),
                                                    rownames(seurat@assays$ISO@data))])
  DefaultAssay(seurat) <- "ISO"
  if (length (x) != 0) {
    p = SpatialFeaturePlot(seurat, features = x, pt.size.factor=1.3, alpha = c(1, 1),ncol = )
    plot_list[[i]] <- p
  }}


filename <- "~/onedrive/Work/phD/phd_project/SiT/results/after_ln_minion/malignant_isoform_RPS_spatial_distribution.pdf"
pdf(filename, width=20, height=10)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()
