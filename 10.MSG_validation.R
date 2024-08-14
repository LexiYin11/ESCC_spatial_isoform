####### 0. Library import #######
library(survival)
library(survminer)
####### 1. Data import #######

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
metadata_ori <- read.csv(file = "rawdata/baiduII_clinical_data.csv", header = TRUE)

metadata <- metadata_ori[,c("sample.ID............BDESCC2..", "Gender",
                               "Location.1","Smoking.history","TNM.stage.the.Eighth.Edition.",
                               "Grade","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.")]
colnames(metadata)[1] <- "Sample_ID"
metadata = metadata[order(metadata$Sample_ID),]
metadata <- metadata[metadata$Sample_ID %in% colnames(count), ]
dim(metadata)


#### mutation data 
#mut_baidu_IandII_combine  <- readRDS( "rawdata/mut_baidu_IandII_combine.rds")

####### 1.2. MSG_gene_list ##############
up_msg_genes <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/function/up_malignant_genes.csv",row.names = 1)
down_msg_genes <- read.csv("~/onedrive/Work/phD/phd_project/SiT/results/function/down_malignant_genes.csv",row.names = 1)


## up selection (archive)
up_genes_ensemble <- HGNC[HGNC$Approved.symbol %in% rownames(up_msg_genes),]
up_genes_ensemble <- up_genes_ensemble$Ensembl.gene.ID

up_msg_genes_count <- count[up_genes_ensemble,]
no_metastasis_sample <- metadata[ metadata$recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free. == "0","Sample_ID"]
no_metastasis <- up_msg_genes_count[, as.character(no_metastasis_sample)]

metastasis_sample <- metadata[ metadata$recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free. == "1","Sample_ID"]
metastasis <- up_msg_genes_count[, as.character(metastasis_sample)]

metastasis_rowmeans_sub <- data.frame(id = rownames(metastasis), sub = rowMeans(metastasis) - rowMeans(no_metastasis),
                                      var_meta = rowRanges(as.matrix(metastasis)),var_no_meta = rowRanges(as.matrix(no_metastasis)) )
metastasis_rowmeans_sub <- metastasis_rowmeans_sub[metastasis_rowmeans_sub$sub > 0,]
metastasis_rowmeans_sub <- metastasis_rowmeans_sub[metastasis_rowmeans_sub$var_meta < 10000,]
metastasis_rowmeans_sub <- metastasis_rowmeans_sub[metastasis_rowmeans_sub$var_no_meta < 1000,]

metastasis_rowmeans_sub <- metastasis_rowmeans_sub[order(metastasis_rowmeans_sub$sub,decreasing = T),]

up_selected_genes <- metastasis_rowmeans_sub$id[1:10]
  #selected_genes[!grepl("NA", selected_genes)]

## down selection

down_genes_ensemble <- HGNC[HGNC$Approved.symbol %in% rownames(down_msg_genes),]
down_genes_ensemble <- down_genes_ensemble$Ensembl.gene.ID
#down_genes_ensemble <- down_genes_ensemble[!grepl("NA.1", down_genes_ensemble)]
down_msg_genes_count <- count[down_genes_ensemble,]

no_metastasis_sample <- metadata[ metadata$recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free. == "0","Sample_ID"]
no_metastasis <- down_msg_genes_count[, as.character(no_metastasis_sample)]

metastasis_sample <- metadata[ metadata$recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free. == "1","Sample_ID"]
metastasis <- down_msg_genes_count[, as.character(metastasis_sample)]

selected_genes <- rownames(metastasis[(rowMeans(metastasis) - rowMeans(no_metastasis)) > 1000 ,])
down_selected_genes <- selected_genes[!grepl("NA", selected_genes)]

msg_genes_list <- list()
msg_genes_list$up_MSG <- up_selected_genes
msg_genes_list$down_MSG <- down_selected_genes
# msg_genes <- c(rownames(up_msg_genes), rownames(down_msg_genes))
# genes_ensemble <- HGNC[HGNC$Approved.symbol %in% msg_genes,]
# genes_ensemble <- genes_ensemble$Ensembl.gene.ID
# msg_genes_list <- list()
# msg_genes_list$MSG <- genes_ensemble
####### 1.3. GSVA  #######
gsva.res <- gsva(as.matrix(count), msg_genes_list, method = "gsva")
dim(gsva.res)
####### 1.4. Clinical information (Fig S9B) ######
gsva.res <- gsva.res[,as.character(metadata$Sample_ID)]
## Metastasis 
metastasis_df <- data.frame(scores = gsva.res["down_MSG",], meta_condition = metadata$recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.)
metastasis_df <- na.omit(metastasis_df )
metastasis_df <- metastasis_df[metastasis_df$meta_condition != "loss",]
pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/new_selected_down_MSG_scores_metastasis_boxplot.pdf",width = 10, height = 8)
ggboxplot(metastasis_df,x = "meta_condition", y = "scores",color = "meta_condition", palette = "jco") + 
  stat_compare_means(comparisons = list(c("0","1"))) 
dev.off()

## Grade 
grade_df <- data.frame(scores = gsva.res["down_MSG",], meta_condition = metadata$Grade)
grade_df <- na.omit(grade_df )
#grade_df$meta_condition <- ifelse(grade_df$meta_condition == "G3", "G3", "G2&G1")
grade_df$meta_condition <- factor(grade_df$meta_condition, levels = c("G1","G2","G3"))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/down_MSG_scores_grade_boxplot.pdf",width = 10, height = 8)

ggboxplot(grade_df,x = "meta_condition", y = "scores",color = "meta_condition", palette = "jco") + 
  stat_compare_means(comparisons = list(c("G1","G3"),c("G1","G2"),c("G2","G3"))) 
dev.off()

## TNM
TNM_df <- data.frame(scores = gsva.res["down_MSG",], meta_condition = metadata$TNM.stage.the.Eighth.Edition.)
TNM_df$meta_condition <- ifelse(TNM_df$meta_condition %in% c("IA",  "IB"), "I",
                     ifelse(TNM_df$meta_condition %in% c("IIA",  "ⅡA","IIB"), "II",
                            ifelse(TNM_df$meta_condition %in% c("IIIA", "IIIB"), "III","IV")))
TNM_df <- na.omit(TNM_df)
TNM_df$meta_condition <- factor(TNM_df$meta_condition, levels = c("I", "II","III", "IV"))
#grade_df <- grade_df[grade_df$meta_condition != "loss",]
pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/down_MSG_scores_TNM_boxplot.pdf",width = 10, height = 8)
ggboxplot(TNM_df,x = "meta_condition", y = "scores",color = "meta_condition", palette = "jco") + 
  stat_compare_means(comparisons = list(c("I","II"),
                                        c("II","III"),
                                        c("III","IV"),
                                        c("I", "III"),
                                        c("I", "IV"))) 
dev.off()


####### 1.5 Survival plot (Fig S9C)######
gsva_res_down <- gsva.res[2,]
summary(gsva_res_down)
patient_high_score <- names(gsva_res_down[gsva_res_down >= 0])
patient_low_score <- names(gsva_res_down[gsva_res_down < 0])
length(patient_high_score)
length(patient_low_score)

metadata_survival<- metadata_ori[,c("sample.ID............BDESCC2..", "overall.survival...............1..dead..0.alive.",
                                "overall.survival........time","Gender","Age", "Location","Family.History.of.ESCC","Family.History.of.other.cancers","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.","TNM.stage.the.Eighth.Edition.","TNM.stage..the.7th.Edition.","Grade","Smoking.history","Drinking.history")]
colnames(metadata_survival)<- c("Sample_ID","Survival_status","Survival_time",  "Gender", "Age", "Location", "ESCC_family_history", "Other_cancer_family_history", "Recurrence_status","TNM_stage_8","TNM_stage_7","Grade","smoking","drinking")
metadata_survival = metadata_survival[order(metadata_survival$Sample_ID),]
metadata_survival <- metadata_survival[metadata_survival$Survival_status != "loss",]
metadata_survival$Survival_status <- as.numeric(metadata_survival$Survival_status)
metadata_survival$Survival_time <- as.numeric(metadata_survival$Survival_time)
metadata_survival$Age <- as.numeric(metadata_survival$Age)
metadata_survival$MSG_score <- ifelse(metadata_survival$Sample_ID %in% patient_high_score,"High", "Low")

km_MSG_score_fit <- survfit(Surv(Survival_time, Survival_status) ~ MSG_score, data=metadata_survival)

pdf("~/onedrive/Work/phD/phd_project/SiT/results/function/down_MSG_scores_survival.pdf",width = 10, height = 6)
ggsurvplot(km_MSG_score_fit, linetype = "strata", break.time.by = 250, 
           palette = c("#E7B800", "#2E9FDF"), conf.int = T, pval = T, data = metadata_survival, risk.table = T)

ggsave("~/onedrive/Work/phD/phd_project/SiT/results/function/down_MSG_scores_survival.pdf")


#### save MSG score gene set 
final_msg_score_gene_df <- HGNC[HGNC$Ensembl.gene.ID %in% down_selected_genes,]

write.csv(final_msg_score_gene_df,"~/onedrive/Work/phD/phd_project/SiT/results/function/final_msg_score_gene_df.csv")






                              