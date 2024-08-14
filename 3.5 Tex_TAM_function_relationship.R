######## 0. Library import #######
library(ggstatsplot)
celltype_cols <- c( "#51c4c2", "#0d8a8c","#4583b3",
                    "#c63596" , "#be86ba", "#8b66b8",
                    "#f78e26",
                    "#f172ad", "#f7afb9",
                    "#4068b2","#512a93","#223271"
)
celltype_cols <- colorRampPalette(celltype_cols)(17)
scales::show_col(alpha(celltype_cols,1))
names(celltype_cols) <- levels(sce_anno$celltype)
######## 1. Data import #######
TAM_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/TAM_only_sce_anno.rds")

tcells_only_sce_anno <- readRDS("~/onedrive/Work/phD/phd_project/SiT/results/sc/merged/remove_doublet/tcells_only_sce_anno.rds")

####  3. Addmodule score #####
cyto_module <- c("GZMA", "GZMB", "GZMH", "GZMM", "GZMK", "GNLY", "PRF1" ,"CTSW")
infl_module <- c("CCL2", "CCL3", "CCL4",
                 "CCL5", "CXCL10", "CXCL9", "IL1B", "IL6", "IL7", "IL15" , "IL18")
hla_dep_inhi_receptor <- c("KIR2DL1", "KIR2DL3",
                           "KIR3DL1", "KIR3DL2", "LILRB1", "LAG3")
hla_indep_inhi_receptor <- c("PDCD1", "SIGLEC7", "CD300A", "CD96", #"IL1RAPL1", 
                             "TIGIT",
                             "HAVCR2")
hla_dep_act_receptor <- c("KIR2DL4", "CD160", "KLRC2")
hla_indep_act_receptor <- c("NCR3", "NCR1", "KLRK1", "CRTAM", "FCGR3A")
## score
tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = cyto_module, name = "Cytotoxicity score")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[16] <- "Cytotoxicity.score"

tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = infl_module, name = "Inflammatory score")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[24] <- "Inflammatory.score"

tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = hla_dep_inhi_receptor, name = "HLA-dependent inhibitory receptors")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[35] <- "HLA.dependent.inhibitory.receptors"

tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = hla_indep_inhi_receptor, name = "HLA-independent inhibitory receptors")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[41] <- "HLA.independent.inhibitory.receptors"

tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = hla_dep_act_receptor, name = "HLA-dependent activating receptors")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[47] <- "HLA.dependent.activating.receptors"

tcells_only_sce_anno <- AddModuleScore(tcells_only_sce_anno, features = hla_indep_act_receptor, name = "HLA-independent activating receptors")
colnames(tcells_only_sce_anno@meta.data)
colnames(tcells_only_sce_anno@meta.data)[50] <- "HLA.independent.activating.receptors"

#### 4. Add module score visualization (Fig 2H,I) #####
## Ridgeplot
pdf("/Users/lexiiiii/onedrive/Work/phD/phd_project/XBP1_innate/results/nkcells/scores_plot.pdf",width = 15)
VlnPlot(tcells_only_sce_anno, features = c("Cytotoxicity.score", "Inflammatory.score",
                                        "HLA-dependent inhibitory receptors",
                                        "HLA-dependent activating receptors",
                                        "HLA-independent activating receptors"),
        group.by = "celltype")
RidgePlot(tcells_only_sce_anno, features = c("Cytotoxicity.score", "Inflammatory.score",
                                          "HLA-dependent inhibitory receptors",
                                          "HLA-dependent activating receptors",
                                          "HLA-independent activating receptors"),
          group.by = "celltype")
VlnPlot(tcells_only_sce_anno, features = c("Cytotoxicity.score", "Inflammatory.score",
                                        "HLA-dependent inhibitory receptors",
                                        "HLA-dependent activating receptors",
                                        "HLA-independent activating receptors"),
        group.by = "celltype")
RidgePlot(tcells_only_sce_anno, features = c("Cytotoxicity.score", "Inflammatory.score",
                                          "HLA-dependent inhibitory receptors",
                                          "HLA-dependent activating receptors",
                                          "HLA-independent activating receptors"),
          group.by =  "celltype")
dev.off()

### dot plot cyto infla
nkcell_metadata <- tcells_only_sce_anno@meta.data
cyto_inf_df <- nkcell_metadata %>% group_by(celltype) %>% summarize(mean_cyto = mean(Cytotoxicity.score), mean_infl = mean(Inflammatory.score)) 
cyto_inf_df <- as.data.frame(cyto_inf_df)

p1 <- ggplot(data = cyto_inf_df ,
             aes(x=mean_cyto, y=mean_infl)) +
  geom_point(aes(fill=celltype), shape=21, size = 15) +
  #geom_smooth(se=FALSE, color=pal[1]) +
  scale_fill_manual(values =celltype_cols) + geom_text(aes(label = celltype), size = 3,colour = "white") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'T cell cluster') + xlab("Cytotoxicity score") + ylab("Inflammatory score")+
  theme_minimal() +
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'),
        legend.position = "none")
### dot plot HLAa
nkcell_metadata <- tcells_only_sce_anno@meta.data
cyto_inf_df <- nkcell_metadata %>% group_by(celltype) %>% summarize(mean_dep = mean(HLA.dependent.activating.receptors), mean_indep = mean(HLA.independent.activating.receptors)) 
cyto_inf_df <- as.data.frame(cyto_inf_df)

p2 <- ggplot(data = cyto_inf_df ,
             aes(x=mean_dep, y=mean_indep)) +
  geom_point(aes(fill=celltype), shape=21, size = 15) +
  #geom_smooth(se=FALSE, color=pal[1]) +
  scale_fill_manual(values = celltype_cols) + geom_text(aes(label = celltype), size = 3,colour = "white") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'T cell cluster') + xlab("HLA-dependent activating receptors score") + ylab("HLA-independent activating receptors score")+
  theme_minimal() + 
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'),
        legend.position = "none")


### dot plot HLAi
nkcell_metadata <- tcells_only_sce_anno@meta.data
cyto_inf_df <- nkcell_metadata %>% group_by(celltype) %>% summarize(mean_dep = mean(HLA.dependent.inhibitory.receptors), mean_indep = mean(HLA.independent.inhibitory.receptors)) 
cyto_inf_df <- as.data.frame(cyto_inf_df)

p3 <- ggplot(data = cyto_inf_df ,
             aes(x=mean_dep, y=mean_indep)) +
  geom_point(aes(fill=celltype), shape=21, size = 15) +
  #geom_smooth(se=FALSE, color=pal[1]) +
  scale_fill_manual(values = celltype_cols) + geom_text(aes(label = celltype), size = 3,colour = "white") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'T cell cluster') + xlab("HLA-dependent inhibitory receptors score") + ylab("HLA-independent inhibitory receptors score")+
  theme_minimal() + 
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'),
        legend.position = "none")
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/Tcell_function_cluster_scores_plot.pdf",width = 5, height = 5)
print(p1)
print(p2)
print(p3)
dev.off()




#### 5. MIF correlation (Fig 2G)######

mif_df <- tcells_only_sce_anno@assays$RNA@counts["MIF",]
mif_df <- mif_df[rownames(tcells_only_sce_anno@meta.data)]
merged_mf_df <- cbind(tcells_only_sce_anno@meta.data, mif_df)

cor(merged_mf_df[,c("Cytotoxicity.score","Inflammatory.score","HLA.independent.activating.receptors", "HLA.dependent.activating.receptors",
                    "HLA.independent.inhibitory.receptors", "HLA.dependent.inhibitory.receptors","mif_df")])
## 没啥correlation

#### 6. MIF expression #####
cyto_mif_df <- merged_mf_df %>% group_by(celltype, orig.ident) %>% summarize(mean_cyto = mean(Cytotoxicity.score), mean_mif = mean(mif_df)) 
cyto_mif_df <- as.data.frame(cyto_mif_df)

ggplot(data = cyto_mif_df ,
             aes(x=mean_cyto, y=mean_mif)) +
  geom_point(aes(fill=celltype), shape=21, size = 15) +
  #geom_smooth(se=FALSE, color=pal[1]) +
  scale_fill_manual(values =celltype_cols) + geom_text(aes(label = celltype), size = 3,colour = "white") +
  #scale_y_continuous(breaks = seq(50000, 200000, 50000)) +
  labs(title = 'T cell cluster') + xlab("Cytotoxicity score") + ylab("MIF score")+
  theme_minimal() +
  theme(plot.title =  element_text(color = "grey20", size = 20, face = 'bold', hjust = 0),
        plot.subtitle =  element_text(color = "grey20", size = 14, hjust = 0),
        plot.caption  =  element_text(color = "grey20", size = 12, face = 'italic', hjust = 1),
        plot.background = element_rect(fill = 'white', color='white'),
        axis.text = element_text(color = "grey20", size = 12),  
        axis.title = element_text(color = "grey20", size = 16, face = 'bold'),
        legend.position = "none")

##### only CD8+ 
merged_mf_df_tex <- merged_mf_df %>% filter(celltype == "CD8+ Effctor T")
cor(merged_mf_df_tex[,c("Cytotoxicity.score","Inflammatory.score","HLA.independent.activating.receptors", "HLA.dependent.activating.receptors",
                    "HLA.independent.inhibitory.receptors", "HLA.dependent.inhibitory.receptors","mif_df")])

## sep high MIF vs Low MIF
merged_mf_df_tex$MIF_label <- ifelse(merged_mf_df_tex$mif_df >= mean(merged_mf_df_tex$mif_df),"High MIF expression", "Low MIF expression")
table(merged_mf_df_tex$MIF_label)

# boxplot
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/MIF_function_Tcell_boxplot.pdf")
for (score_type in colnames(merged_mf_df_tex)[1:6]){
  cur_df <- merged_mf_df_tex[,c(score_type, "MIF_label")]
  colnames(cur_df)[1] <- c("function_score")
  p1 <- ggplot(cur_df, aes( MIF_label, function_score, fill=MIF_label))+
    geom_boxplot(outlier.shape=NA,fatten = 2)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width=0.2, size=1, color='midnightblue', alpha=0.8)+
    ylab(score_type)+xlab('') +
    #scale_fill_manual(values = c("Female" ="#CD534CFF","Male" = "#7AA6DCFF")) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          strip.text = element_text(size = 12),axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          legend.position = "none", axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = 'bold'), 
          plot.title = element_text(hjust =0.5)) + stat_compare_means(method = "t.test")
  print(p1)
}
dev.off()

###### Tex MIF cells vs function #####
merged_mf_df_tex <- merged_mf_df_tex[,c("Cytotoxicity.score","Inflammatory.score","HLA.independent.activating.receptors", "HLA.dependent.activating.receptors",
                                        "HLA.independent.inhibitory.receptors", "HLA.dependent.inhibitory.receptors","mif_df","MIF_label")]


mat <- merged_mf_df_tex %>% group_by(MIF_label) %>% summarize(Cytotoxicity.score = mean(Cytotoxicity.score),
                                                              Inflammatory.score = mean(Inflammatory.score),
                                                              HLA.independent.activating.receptors = mean(HLA.independent.activating.receptors),
                                                              HLA.dependent.activating.receptors = mean(HLA.dependent.activating.receptors),
                                                              HLA.independent.inhibitory.receptors = mean(HLA.independent.inhibitory.receptors),
                                                              HLA.dependent.inhibitory.receptors = mean(HLA.dependent.inhibitory.receptors),
                                                              MIF = mean(mif_df))
mat <- as.matrix(mat[,2:6])
mat <- as.matrix(t(mat))
colnames(mat) <- c("High MIF expression", "Low MIF expression")
mat <- as.matrix(mat)
mat <- scale_mat(mat, "column")

dim(mat)
col_fun <- colorRamp2(c(1, 0, -1), c("#fc5d59", "#ffffbf","#91bfdb" ))
pdf("~/onedrive/Work/phD/phd_project/SiT/results/sc/res_1/M1_function_gene_TAM_circos.pdf")
circos.par(gap.after=10)
circos.heatmap(as.matrix(mat),  col = col_fun,rownames.side = "outside")
               cluster = T,track.height =0.4, rownames.cex = 1.2,
               cell_width = 10) 
circos.clear()
dev.off()
