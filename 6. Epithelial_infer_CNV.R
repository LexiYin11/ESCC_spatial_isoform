library(infercnv)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
library(ggsci)
library(SpatialInferCNV)
######## 1. extract epithelial cells only #####
#要求不大 两个epithlial group 有区别就好
after <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after/p1_after.seurat.rds")
before <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/before/Result_X101SC21081081-Z01-J001/3.Seurat/P1_before/P1_before.seurat.rds")
after_ln <- readRDS("~/onedrive/Work/phD/phd_project/SiT/rawData/after/Result_X101SC21081081-Z01-J002\ 2/3.Seurat/p1_after_LN/p1_after_LN.seurat.rds")



## change current data
spdata <- after_ln

######## 2. matrix and file prep  ######
dat = as.data.frame(GetAssayData(spdata, slot='counts' ,assay = 'Spatial'))#表达矩阵
dat[1:4,1:4]
dim(dat)

groupinfo = data.frame(v1=colnames(dat),v2=spdata$celltype)#样本注释文件
head(groupinfo)

geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
head(geneInfor)

dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)

## save
setwd("~/onedrive/Work/phD/phd_project/SiT/results/after_ln/inferCNV/")
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

######### 3. Run inferCNV ######
expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt'
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=unique(after_ln$celltype)[unique(after_ln$celltype) != "Malignant Epithelial cells, Tcells, Macrophage"])  ## 不知道哪个是正常的

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir=  './' ,  #选择输入文件夹
                              cluster_by_groups=T,   # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                              plot_steps=F,num_threads = 30 ,
                              denoise=T, HMM=F,output_format = "pdf")


######### 4. Results Analysis (Fig S8A) ######
## after
setwd("~/onedrive/Work/phD/phd_project/SiT/results/after/inferCNV/SCT/ref_groups_NULL")
### get groupings
group <- read.table("infercnv.observation_groupings.txt", sep = " ", header = T)
### calculate score
obs_ori <- read.table("infercnv.observations.txt", header = T, check.names = F)

max(obs_ori)
min(obs_ori)
obs <- obs_ori

obs[obs > 0.6 & obs < 0.9] <- 2   #把0.6-0.7定义为数值2  后面依此类推
obs[obs >= 0.9 & obs < 0.95] <- 1
obs[obs >= 0.95 & obs < 0.05] <- 0
obs[obs >= 0.05 & obs <= 1] <- 1
obs[obs > 1 & obs <= 1.2] <- 2

#table(obs)
scores <- as.data.frame(colSums(obs))
scores$cluster <- group$Annotation.Group
colnames(scores) <- c("scores", "cluster")
summary(scores)
ggplot(scores, aes(x=cluster, y=scores, fill=factor(cluster)))+
  geom_boxplot(outlier.shape=NA)+
  #scale_fill_simpsons() +
  labs(fill = "Cluster")



  