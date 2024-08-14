library(survival)
library(survminer)
library(EnhancedVolcano)
library(enrichplot)
library(clusterProfiler)
######### 0. Markers #######
#### Between C1QC+ TAM and FCN1+ TAM
#TAM_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/res_0.5_TAM_only_sce_anno_anno.rds")
TAM_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_only_sce_anno.rds")

c1qc_TAM_markers <- FindMarkers(TAM_only_sce_anno, group.by = "tam_label",ident.1 = "C1QC+ TAM", ident.2 = "FCN1+ TAM")
dim(c1qc_TAM_markers)
#write.csv(c1qc_TAM_markers, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_TAM_markers.csv", row.names = F)
write.csv(c1qc_TAM_markers, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_TAM_markers.csv", row.names = F)

c1qc_TAM_markers$gene <- rownames(c1qc_TAM_markers)
c1qc_TAM_markers = c1qc_TAM_markers %>% select(gene, everything()) %>% subset(p_val<0.05)
c1qc_TAM_markers <- c1qc_TAM_markers[c1qc_TAM_markers$p_val <= 0.05,]
#先把基因这一列放在第一列，然后选取p值小于0.05的行(结果行数不变，说明都挺好的)
c1qc_TAM_up = c1qc_TAM_markers[c1qc_TAM_markers$avg_log2FC >=1,]
  c1qc_TAM_markers %>% select(gene, everything()) %>% subset(avg_log2FC >= 1)
c1qc_TAM_down = c1qc_TAM_markers%>% select(gene, everything()) %>% subset(avg_log2FC <= -1)

write.csv(c1qc_TAM_up, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_TAM_up.csv", row.names = F)
write.csv(c1qc_TAM_down, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_TAM_down.csv", row.names = F)

######### 1. Volcano plot （Fig S5B) #######
#### Between C1QC+ TAM and FCN1+ TAM
keyvals <- rep('black', nrow(c1qc_TAM_markers))

# set the base name/label as 'Mid'
names(keyvals) <- rep('Mid', nrow(c1qc_TAM_markers))

# fold change > 1.5 & p-value < 0.0001 为高表达
keyvals[which(c1qc_TAM_markers$avg_log2FC > 1 & c1qc_TAM_markers$p_val_adj<0.001)] <- 'gold'
names(keyvals)[which(c1qc_TAM_markers$avg_log2FC > 1 & c1qc_TAM_markers$p_val_adj<0.001)] <- 'high'

# fold change < -1.5 & p-value < 0.0001为低表达
keyvals[which(c1qc_TAM_markers$avg_log2FC < -1 & c1qc_TAM_markers$p_val_adj<0.001)] <- 'royalblue'
names(keyvals)[which(c1qc_TAM_markers$avg_log2FC < -1 & c1qc_TAM_markers$p_val_adj<0.001)] <- 'low'
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/c1qc_diff_markers_EnhancedVolcan.pdf")
EnhancedVolcano(c1qc_TAM_markers,
                lab = rownames(c1qc_TAM_markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-4, 4),
                pCutoff = 0.001,
                FCcutoff = 1,
                labSize = 3,
                pointSize = 3,
                captionLabSize = 8,
                titleLabSize = 8,
                subtitleLabSize = 6,
                axisLabSize = 8,
                legendPosition = 'right',
                colCustom = keyvals)
                
dev.off()

######### 2. Baidu II data #######

### count
count <- read.table("rawdata/baiduII_gene_count_.txt",header = T, row.names = 1)
nrow(count)
nrow(na.omit(count))
row_list <- c()

for (rowname in rownames(count)) {
  row_list <- c(row_list , strsplit(rowname, "_")[[1]][1])
}
rownames(count) <- row_list
colnames(count) <- gsub("X","",colnames(count))
HGNC <- read.delim("rawdata/HGNC_results.txt")
### only keep tumor
count <- count[,grepl(".T",colnames(count))]
colnames(count) <- gsub("T","",colnames(count))
colnames(count) <- as.integer(colnames(count))
## clinical data 
metadata_baiduII_ori <- read.csv(file = "rawdata/baiduII_clinical_data.csv", header = TRUE)

metadata_BaiduII <- metadata_baiduII_ori[,c("sample.ID............BDESCC2..", "Gender",
                                            "Smoking.history",  "Grade","TNM.stage.the.Eighth.Edition.","Age", "T","N","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.")]
colnames(metadata_BaiduII) <- c("Tumor_Sample_Barcode","Gender","Smoking status","Grade","TNM stage", "Age","T", "N","Metastasis status")
metadata_BaiduII$`TNM stage` <- ifelse(metadata_BaiduII$`TNM stage` %in% c("IIIA","IIIB"),3,
                                       ifelse(metadata_BaiduII$`TNM stage` %in% c("IIA","ⅡA","IIB"),2,
                                              ifelse(metadata_BaiduII$`TNM stage` == "IB",1,
                                                     ifelse(metadata_BaiduII$`TNM stage` == "ⅣA",4,"stop"))))


#### mutation data 
common_colnames <- intersect(colnames(mut_baiduI),colnames(mut_baiduII))
pre_mut_baidu_IandII_combine <- rbind(mut_baiduI[,common_colnames] ,
                                      mut_baiduII[,common_colnames] )
write.csv(pre_mut_baidu_IandII_combine, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/pre_mut_baidu_IandII_combine.csv")
pre_mut_baidu_IandII_combine <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/pre_mut_baidu_IandII_combine.csv")

####### 2.1 GSVA #######
c1qc_TAM_up_en <- HGNC[HGNC$Approved.symbol %in% c1qc_TAM_up$gene,]
c1qc_TAM_up_en <- c1qc_TAM_up_en$Ensembl.gene.ID
c1qc_TAM_up_en <- list(c1qc_TAM_up = c1qc_TAM_up_en)
gsva.c1qc_TAM_markers <- gsva(as.matrix(count), c1qc_TAM_up_en, method = "gsva")
dim(gsva.c1qc_TAM_markers)
gsva.c1qc_TAM_markers <- as.data.frame(t(gsva.c1qc_TAM_markers))
gsva.c1qc_TAM_markers$sample <- rownames(gsva.c1qc_TAM_markers)
summary(gsva.c1qc_TAM_markers)

patient_high_score <- gsva.c1qc_TAM_markers[gsva.c1qc_TAM_markers$c1qc_TAM_up >= mean(gsva.c1qc_TAM_markers$c1qc_TAM_up ),"sample"]
patient_low_score <- gsva.c1qc_TAM_markers[gsva.c1qc_TAM_markers$c1qc_TAM_up < mean(gsva.c1qc_TAM_markers$c1qc_TAM_up ),"sample"]
length(patient_high_score)
length(patient_low_score)


####### 2.2 exp numerical clinical information ######
gsva.c1qc_TAM_markers <- gsva.c1qc_TAM_markers[as.character(metadata_BaiduII$Tumor_Sample_Barcode),]
clinical_gsva_df <- merge(metadata_BaiduII,gsva.c1qc_TAM_markers, by.x = "Tumor_Sample_Barcode", by.y = "sample")
clinical_gsva_df <- clinical_gsva_df[clinical_gsva_df$`TNM stage` !="stop",]
clinical_gsva_df <- clinical_gsva_df[clinical_gsva_df$`Metastasis status` != "loss",]
clinical_gsva_df$`TNM stage` <- factor(clinical_gsva_df$`TNM stage`, levels = c("1", "2", "3","4"))
clinical_gsva_df$`Smoking status` <- factor(clinical_gsva_df$`Smoking status`, levels = c("never","light", "moderate", "heavy"))
clinical_gsva_df$Age_group <- ifelse(clinical_gsva_df$Age <= 50, "<= 50"," >50")
        #                             ifelse(clinical_gsva_df$Age <= 65, "> 40 & <= 70",">70"))
clinical_gsva_df$N_group <- ifelse(clinical_gsva_df$N %in% c("1","2"), "1 & 2", "3 & 4")

pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/clinical_info/numerical_c1qc_diff_logFC1_clinical_info.pdf")
ggboxplot(clinical_gsva_df,x = "Metastasis status", y = "c1qc_TAM_up",color = "Metastasis status", palette = "jco") +
  stat_compare_means(comparisons = list(c("0","1")))


ggboxplot(clinical_gsva_df,x = "Grade", y = "c1qc_TAM_up",color = "Grade", palette = "jco") +
  stat_compare_means(comparisons = list(c("G1","G3"),c("G1","G2"),c("G2","G3")))

ggboxplot(clinical_gsva_df,x = "TNM stage", y = "c1qc_TAM_up",color = "TNM stage", palette = "jco") +
  stat_compare_means(comparisons = list(c("1","2"),
                                        c("2","3"),
                                        c("3","4"),
                                        c("1", "3"),
                                        c("1", "4")))
ggboxplot(clinical_gsva_df,x = "Smoking status", y = "c1qc_TAM_up",color = "Smoking status", palette = "jco") +
  stat_compare_means(comparisons = list(c("never","light"),c("never","moderate"),
                                        c("never","heavy")))
ggboxplot(clinical_gsva_df,x = "Gender", y = "c1qc_TAM_up",color = "Gender", palette = "jco") +
  stat_compare_means()

ggboxplot(clinical_gsva_df,x = "T", y = "c1qc_TAM_up",color = "T", palette = "jco") +
  stat_compare_means(comparisons = list(c("1","2"), c("1", "3")))

ggboxplot(clinical_gsva_df,x = "N", y = "c1qc_TAM_up",color = "N", palette = "jco") +
  stat_compare_means(comparisons = list(c("1","2"), c("1", "3")))

ggboxplot(clinical_gsva_df,x = "N_group", y = "c1qc_TAM_up",color = "N_group", palette = "jco") +
  stat_compare_means()

ggboxplot(clinical_gsva_df,x = "Age_group", y = "c1qc_TAM_up",color = "Age_group", palette = "jco") +
  stat_compare_means()

 dev.off()

 
####### 2.3 exp categorical clinical information  #########
metadata_BaiduII$group <- ifelse(metadata_BaiduII$Tumor_Sample_Barcode %in% patient_high_score,"High", "Low")
table(metadata_BaiduII$group)
comparison_all<- function(metadata, file_name){
  ### Grade
  p1_grade = ggbarstats(metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata
  new_metadata$Grade <- ifelse(metadata$Grade %in% c("G1"), "G1", "Others")
  p2_grade = ggbarstats(new_metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$Grade <- ifelse(metadata$Grade %in% c("G3"), "G3", "Others")
  p3_grade = ggbarstats(new_metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### TNM
  p1_TNM = ggbarstats(metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`TNM stage` <- ifelse(metadata$`TNM stage` %in% c("4"), "4", "Others")
  p2_TNM = ggbarstats(new_metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`TNM stage` <- ifelse(metadata$`TNM stage` %in% c("4","3"), "4,3", "Others")
  p3_TNM = ggbarstats(new_metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### T
  p1_T = ggbarstats(metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`T` <- ifelse(metadata$`T` %in% c("4"), "4", "Others")
  p2_T = ggbarstats(new_metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`T` <- ifelse(metadata$`T` %in% c("4","3"), "4,3", "Others")
  p3_T = ggbarstats(new_metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### N
  p1_N = ggbarstats(metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`N` <- ifelse(metadata$`N` %in% c("3"), "3", "Others")
  p2_N = ggbarstats(new_metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`N` <- ifelse(metadata$`N` %in% c("2","3"), "3,2", "Others")
  p3_N = ggbarstats(new_metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### Age
  new_metadata <- metadata
  new_metadata$Age <- ifelse(metadata$`Age` >= 60 , ">= 60", "Others")
  p1_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$Age <- ifelse(metadata$`Age` >= 65, ">= 65", "Others")
  p2_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata  
  new_metadata$Age <- ifelse(metadata$Age <= 60, "<= 60", 
                             ifelse(metadata$Age <= 65, "> 60 & <= 65", 
                                    ifelse(metadata$Age <= 70, "> 65 & <= 70", 
                                           ifelse(metadata$Age <= 80, "> 70 & <= 80","> 80"))))
  p3_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### Metastasis
  p1_Metastasis = ggbarstats(metadata, "Metastasis status", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  
  pdf(paste0("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/clinical_info/",file_name,"baiduI_II_clinical_info_mutation.pdf"),width = 20)
  ggarrange(p1_grade,p2_grade,p3_grade,p1_TNM,p2_TNM,p3_TNM, p1_T,p2_T,p3_T, 
            p1_N,p2_N,p3_N, p1_Age,p2_Age,p3_Age,p1_Metastasis ,ncol = 3)
  dev.off()
}
metadata <- metadata_BaiduII
file_name <- "exp_c1qc_diff_logFC_1_"



####### 2.4 exp Survival plot ######
metadata_survival<- metadata_ori[,c("sample.ID............BDESCC2..", "overall.survival...............1..dead..0.alive.",
                                    "overall.survival........time","Gender","Age", "Location","Family.History.of.ESCC","Family.History.of.other.cancers","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.","TNM.stage.the.Eighth.Edition.","TNM.stage..the.7th.Edition.","Grade","Smoking.history","Drinking.history")]
colnames(metadata_survival)<- c("Sample_ID","Survival_status","Survival_time",  "Gender", "Age", "Location", "ESCC_family_history", "Other_cancer_family_history", "Recurrence_status","TNM_stage_8","TNM_stage_7","Grade","smoking","drinking")
metadata_survival = metadata_survival[order(metadata_survival$Sample_ID),]
metadata_survival <- metadata_survival[metadata_survival$Survival_status != "loss",]
metadata_survival$Survival_status <- as.numeric(metadata_survival$Survival_status)
metadata_survival$Survival_time <- as.numeric(metadata_survival$Survival_time)
metadata_survival$Age <- as.numeric(metadata_survival$Age)
metadata_survival$group <- ifelse(metadata_survival$Sample_ID %in% patient_high_score,"High", "Low")

km_MSG_score_fit <- survfit(Surv(Survival_time, Survival_status) ~ group, data=metadata_survival)
ggsurvplot(km_MSG_score_fit, linetype = "strata", break.time.by = 250, 
           palette = c("#E7B800", "#2E9FDF"), conf.int = T, pval = T, data = metadata_survival, risk.table = F)

ggsave("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/clinical_info/exp_c1qc_diff_logfc_1_survival.pdf")



####### 2.5 mut categorical 临床信息 （baidu I and baidu II) #####
altered_patients <- pre_mut_baidu_IandII_combine[pre_mut_baidu_IandII_combine$Hugo_Symbol %in% c1qc_TAM_up$gene,]
altered_patients  <- unique(altered_patients$Tumor_Sample_Barcode)
length(altered_patients)

metadata_baiduI_BaiduII$group <- ifelse(metadata_baiduI_BaiduII$Tumor_Sample_Barcode %in% altered_patients, "mutated", "unmutated")
table(metadata_baiduI_BaiduII$group)
comparison_all<- function(metadata, file_name){
  ### Grade
  p1_grade = ggbarstats(metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata
  new_metadata$Grade <- ifelse(metadata$Grade %in% c("G1"), "G1", "Others")
  p2_grade = ggbarstats(new_metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$Grade <- ifelse(metadata$Grade %in% c("G3"), "G3", "Others")
  p3_grade = ggbarstats(new_metadata, Grade, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### TNM
  p1_TNM = ggbarstats(metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`TNM stage` <- ifelse(metadata$`TNM stage` %in% c("4"), "4", "Others")
  p2_TNM = ggbarstats(new_metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`TNM stage` <- ifelse(metadata$`TNM stage` %in% c("4","3"), "4,3", "Others")
  p3_TNM = ggbarstats(new_metadata, `TNM stage`, group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### T
  p1_T = ggbarstats(metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata     
  new_metadata$`T` <- ifelse(metadata$`T` %in% c("4"), "4", "Others")
  p2_T = ggbarstats(new_metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`T` <- ifelse(metadata$`T` %in% c("4","3"), "4,3", "Others")
  p3_T = ggbarstats(new_metadata, "T", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### N
  p1_N = ggbarstats(metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`N` <- ifelse(metadata$`N` %in% c("3"), "3", "Others")
  p2_N = ggbarstats(new_metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$`N` <- ifelse(metadata$`N` %in% c("2","3"), "3,2", "Others")
  p3_N = ggbarstats(new_metadata, "N", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### Age
  new_metadata <- metadata
  new_metadata$Age <- ifelse(metadata$`Age` >= 60 , ">= 60", "Others")
  p1_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata          
  new_metadata$Age <- ifelse(metadata$`Age` >= 65, ">= 65", "Others")
  p2_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  new_metadata <- metadata  
  new_metadata$Age <- ifelse(metadata$Age <= 60, "<= 60", 
                             ifelse(metadata$Age <= 65, "> 60 & <= 65", 
                                    ifelse(metadata$Age <= 70, "> 65 & <= 70", 
                                           ifelse(metadata$Age <= 80, "> 70 & <= 80","> 80"))))
  p3_Age = ggbarstats(new_metadata, "Age", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  ### Metastasis
  p1_Metastasis = ggbarstats(metadata, "Metastasis status", group, palette = 'Set2',ggstatsplot.layer = FALSE)
  
  
  pdf(paste0("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/clinical_info/",file_name,"baiduI_II_clinical_info_mutation.pdf"),width = 20)
  ggarrange(p1_grade,p2_grade,p3_grade,p1_TNM,p2_TNM,p3_TNM, p1_T,p2_T,p3_T, 
            p1_N,p2_N,p3_N, p1_Age,p2_Age,p3_Age,p1_Metastasis ,ncol = 3)
  dev.off()
}

metadata <- metadata_baiduI_BaiduII
file_name <- "mut_c1qc_diff_logFC_1_categorical_"



####### 2.6 survival plot for mutation ######
##### mutation survival plot
altered_patients <- pre_mut_baidu_IandII_combine[pre_mut_baidu_IandII_combine$Hugo_Symbol %in% c1qc_TAM_up$gene,]
altered_patients  <- unique(altered_patients$Tumor_Sample_Barcode)
length(altered_patients)


##### baiduI and II metadata
metadata_baiduI_ori <- read.csv(file = "~/onedrive/Work/phD/phd_project/TME_gender/rawdata/ESCC_pheno_508_20180911.csv", header = TRUE)
metadata_baiduI_survival <- metadata_baiduI_ori[,c("tumor",  "Survival.status", "Survival.after.diagnosis.months")]
colnames(metadata_baiduI_survival)<- c("Sample_ID","Survival_status","Survival_time")

metadata_baiduII_ori <- read.csv(file = "~/onedrive/Work/phD/phd_project/TME_gender/rawdata/baiduII_clinical_data.csv", header = TRUE)
metadata_baiduII_survival <- metadata_baiduII_ori[,c("sample.ID............BDESCC2..", "overall.survival...............1..dead..0.alive.",
                                                     "overall.survival........time")]
colnames(metadata_baiduII_survival)<- c("Sample_ID","Survival_status","Survival_time")
metadata_baiduI_baiduII_survival <- rbind(metadata_baiduI_survival, metadata_baiduII_survival)



## filter
metadata_baiduI_survival <- metadata_baiduI_survival[!metadata_baiduI_survival$Survival_status %in% c("#VALUE!","Losstofollow-up","loss"),]
metadata_baiduI_survival$Survival_status <- as.numeric(metadata_baiduI_survival$Survival_status)
metadata_baiduI_survival$Survival_time <- as.numeric(metadata_baiduI_survival$Survival_time)
metadata_baiduI_survival$c1qc_diff_mutated <- ifelse(metadata_baiduI_survival$Sample_ID %in% altered_patients, "c1qc_diff_mutated", "c1qc_diff_unmutated")

metadata_baiduII_survival <- metadata_baiduII_survival[!metadata_baiduII_survival$Survival_status %in% c("#VALUE!","Losstofollow-up","loss"),]
metadata_baiduII_survival$Survival_status <- as.numeric(metadata_baiduII_survival$Survival_status)
metadata_baiduII_survival$Survival_time <- as.numeric(metadata_baiduII_survival$Survival_time)
metadata_baiduII_survival$c1qc_diff_mutated <- ifelse(metadata_baiduII_survival$Sample_ID %in% altered_patients, "c1qc_diff_mutated", "c1qc_diff_unmutated")

metadata_baiduI_baiduII_survival <- metadata_baiduI_baiduII_survival[!metadata_baiduI_baiduII_survival$Survival_status %in% c("#VALUE!","Losstofollow-up","loss"),]
metadata_baiduI_baiduII_survival$Survival_status <- as.numeric(metadata_baiduI_baiduII_survival$Survival_status)
metadata_baiduI_baiduII_survival$Survival_time <- as.numeric(metadata_baiduI_baiduII_survival$Survival_time)
metadata_baiduI_baiduII_survival$c1qc_diff_mutated <- ifelse(metadata_baiduI_baiduII_survival$Sample_ID %in% altered_patients, "c1qc_diff_mutated", "c1qc_diff_unmutated")

table(metadata_baiduI_baiduII_survival$c1qc_diff_mutated)


### survival plot
km_fit_mag <- survfit(Surv(Survival_time, Survival_status) ~ c1qc_diff_mutated, data= metadata_baiduI_baiduII_survival)
summary(km_fit_mag)
ggsurvplot(km_fit_mag, data = metadata_baiduI_baiduII_survival,conf.int = T, pval = T)
ggsave("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/mut_c1qc_diff_logFC_1_ggsurvplot_baiduI_baiduII.pdf")


####### 3. TCGA data 90 cases ######

count <- as.data.frame(read.csv("~/onedrive/Work/phD/phd_project/SiT/rawData/TCGA_counts.csv"))
count <- count[!is.na(count$X),]
count <- count[!duplicated(count$X),]

count %>% distinct(X,.keep_all = T) %>% as.tibble() %>% column_to_rownames(var = "X")
rownames(count) <- count$X
count <- count[,colnames(count)[2:186]]
colnames(count) <- gsub("\\.", "-",colnames(count))
metadata <- as.data.frame(read_csv("~/onedrive/Work/phD/phd_project/SiT/rawData/TCGA_clinical_info.csv"))
count <- count[,metadata$`Sample ID`]
dim(count)
## 90 cases
####### 3.1 GSVA ##### 
gsva.res <- gsva(as.matrix(count), c1qc_TAM_up_en, method = "gsva")
gsva.res  <- as.data.frame(t(gsva.res ))
summary(gsva.res )
gsva.res$Sample_ID <- rownames(gsva.res)
patient_high_score <- gsva.res[gsva.res$c1qc_TAM_up >= mean(gsva.res$c1qc_TAM_up ),"Sample_ID"]
patient_low_score <- gsva.res[gsva.res$c1qc_TAM_up < mean(gsva.res$c1qc_TAM_up ),"Sample_ID"]
length(patient_high_score)
length(patient_low_score)

####### 3.2 exp Survival plot 
metadata_survival <-  metadata[,c("case_submitter_id","days_to_death","vital_status")]
colnames(metadata_survival)<- c("Sample_ID","Survival_time", "Survival_status")

metadata_survival$Survival_status <- ifelse(metadata_survival$Survival_status =="Dead",1,0)
metadata_survival$Survival_time <- ifelse(metadata_survival$Survival_time == "'--","1500",metadata_survival$Survival_time)
metadata_survival$Survival_time <- as.numeric(metadata_survival$Survival_time)
metadata_survival <- merge(metadata_survival, gsva.res,by="Sample_ID" )
metadata_survival$group <- ifelse(metadata_survival$Sample_ID %in% patient_high_score,"High", "Low")
km_tcga_fit <- survfit(Surv(Survival_time, Survival_status) ~ group, data=metadata_survival)

ggsurvplot(km_tcga_MSG_score_fit, linetype = "strata", break.time.by = 250, 
           palette = c("#E7B800", "#2E9FDF"), conf.int = T, pval = T, data = metadata_survival, risk.table = F)
ggsave("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/clinical_info/exp_TCGA_c1qc_diff_logfc_2_survival.pdf")



##### 5. GSVA of TAM function (celltype and orig.ident) #####
## gene list
tam_function_genesets <- list(M1 = m1m2_marker$M1, M2 = m1m2_marker$M2,
                              phagocytosis = m1m2_marker$Phagocytosis,
                              Angio = m1m2_marker$Angiogenesis)

modify_tree(tam_function_genesets, leaf = na.omit)
tam_function_genesets <- rapply(tam_function_genesets, na.omit, how = "replace")
tam_function_genesets <- as.list(tam_function_genesets)

tam_function_genesets <- list(M1 = c(tam_function_genesets$M1), M2 = c(tam_function_genesets$M2),
                              phagocytosis = c(tam_function_genesets$phagocytosis),
                              Angio = c(tam_function_genesets$Angio))

expr_mono <- AverageExpression(TAM_only_sce_anno, group.by = c("tam_label", "orig.ident"))
expr_mono <- as.matrix(expr_mono$SCT)
gsva.res_mono <- gsva(expr_mono, tam_function_genesets, method = "gsva")
gsva.res_mono <- as.data.frame(t(gsva.res_mono))
write.csv(gsva.res_mono, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/gsva.res_TAM_celltype_")
### dotplot m1 m2 
gsva_tam_m1_m2_df <- gsva.res_mono
gsva_tam_m1_m2_df$celltype <- sapply(strsplit(rownames(gsva_tam_m1_m2_df),"_"),"[[",1)
gsva_tam_m1_m2_df$patient <- rep(c("after","after_ln","before"),2)
# gsva_tam_m1_m2_df <- data.frame(score = c(gsva.res_mono$M1,gsva.res_mono$M2,
#                                     gsva.res_mono$phagocytosis,gsva.res_mono$Angio),
#                           celltype = c(rep(rownames(gsva.res_mono),4)),
#                           pathway = c(rep("M1",6),rep("M2",6),rep("Phagocytosis",6),rep("Angiogenesis",6)))
p1 <- ggplot(data = gsva_tam_m1_m2_df ,
             aes(x=M1, y=M2)) +
  geom_point(aes(color=patient, shape = celltype), size = 15) +
  scale_color_manual(values = c("#512a93", "#8b66b8","#ffffbf")) +
  scale_shape_manual(values = c(15,16)) +
  geom_text(aes(label = celltype), size = 3,colour = "black") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'TIM') + xlab("M1 score") + ylab("M2 score")+
  theme_minimal() +
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'))

### dotplot phago vs angio 

p2 <- ggplot(data = gsva_tam_m1_m2_df ,
             aes(x=phagocytosis, y=Angio)) +
  geom_point(aes(color=patient, shape = celltype), size = 15) +
  scale_color_manual(values = c("#512a93", "#8b66b8","#ffffbf")) +
  scale_shape_manual(values = c(15,16)) +
  geom_text(aes(label = celltype), size = 3,colour = "black") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'TIM') + xlab("Phagocytosis score") + ylab("Angiogenesis score")+
  theme_minimal() +
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'))

#### GSVA boxplot (all cells) (archive) #######
expr_mono <- TAM_only_sce_anno@assays$SCT@counts
expr_mono <- as.matrix(expr_mono)
gsva.res_mono <- gsva(expr_mono, tam_function_genesets, method = "gsva")
gsva.res_mono <- as.data.frame(t(gsva.res_mono))

#gsva.res_mono <- as.data.frame(t(gsva.res_mono))
write.csv(gsva.res_mono, "~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/gsva.res_TAM_all_cells.csv")

ggboxplot(gsva.res_mono,x = "Grade", y = "c1qc_TAM_up",color = "Grade", palette = "jco") +
  stat_compare_means(comparisons = list(c("G1","G3"),c("G1","G2"),c("G2","G3")))


### 6. add modulescore #####
all_genes <- rownames(TAM_only_sce_anno@assays$SCT@counts)
## score
TAM_only_sce_anno <- AddModuleScore(TAM_only_sce_anno, features = intersect(tam_function_genesets$M1,all_genes), name = "M1 score")
colnames(TAM_only_sce_anno@meta.data)
colnames(TAM_only_sce_anno@meta.data)[21] <- "M1.score"

TAM_only_sce_anno <- AddModuleScore(TAM_only_sce_anno, features = intersect(tam_function_genesets$M2,all_genes), name = "M2 score")
colnames(TAM_only_sce_anno@meta.data)
colnames(TAM_only_sce_anno@meta.data)[63] <- "M2.score"

TAM_only_sce_anno <- AddModuleScore(TAM_only_sce_anno, features = intersect(tam_function_genesets$phagocytosis,all_genes), name = "phagocytosis")
colnames(TAM_only_sce_anno@meta.data)
colnames(TAM_only_sce_anno@meta.data)[90] <- "phagocytosis.score"

TAM_only_sce_anno <- AddModuleScore(TAM_only_sce_anno, features = intersect(tam_function_genesets$Angio,all_genes), name = "Angio")
colnames(TAM_only_sce_anno@meta.data)
colnames(TAM_only_sce_anno@meta.data)[94] <- "Angiogenesis.score"

######
VlnPlot(TAM_only_sce_anno, features = c("M1.score","M2.score","phagocytosis.score",
                                        "Angiogenesis.score"), group.by = "tam_label")


######## 定义 FCN1+ TAM and C1QC+ TAM #####
mono_macro_markers <- c("FCN1", "S100A9","S100A8","FCGR3A", "LST1", "LILRB2", "INHBA",
                        "IL1RN", "CCL4", "NLRP3","EREG", "IL1B","LYVE1",
                        "PLTP","SEPP1","C1QC","C1QB","C1QA")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_definition_marker_dotplot.pdf",width = 13,height = 3)
DotPlot(TAM_only_sce_anno, features = mono_macro_markers , group.by = "tam_label", dot.scale = 12,
        cols = "Spectral")
dev.off()


