library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(seurat_clusters)
library(org.Hs.eg.db)
library(ggpubr)
library(rstatix)

pbmc <- readRDS("result/BMDM.rds")
DefaultAssay(pbmc) <- "RNA"
cluster_color <- readRDS(file = "color/ClusterColor.rds")
col.ls <- cluster_color
names(col.ls) <- names(table(pbmc$seurat_clusters))

NES = list(c('S100A10','FOSL2','SPP1','CAV1','ANXA1','VIM','CD44','SERPINH1',
             'LGALS3','CEBPB','ATF5','LGALS1'))

pbmc = Seurat::AddModuleScore(pbmc, features = NES, name='NES')
p1 <- FeaturePlot(pbmc, features = 'NES1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")
df <- pbmc@meta.data
stat.test <- df %>% t_test(NES1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
df$seurat_clusters <- factor(df$seurat_clusters, level=c("AC","MES","NPC","OPC"))
p2 <- ggboxplot(df, x = "seurat_clusters", y = "NES1", fill = "seurat_clusters", 
                palette = col.ls) + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)

Apoptosis = list(c('ACVR1',	'ACVR1B',	'AKT1',	'AKT2',	'AKT3',	'ATM',	'ATR',
                   'BCL2',	'BCL2L1',	'BMPR1A',	'BMPR1B',	'CAPN1',	'CAPN10',	'CAPN11',	'CAPN12',
                   'CAPN13',	'CAPN14',	'CAPN2',	'CAPN3',	'CAPN5',	'CAPN6',	'CAPN7',	'CAPN8',
                   'CAPN9',	'CASP10',	'CASP12',	'CASP3',	'CASP6',	'CASP7',	'CASP8',	'CASP9',	'CHEK1',
                   'CHEK2',	'CHUK',	'EGFR',	'ERBB2',	'ERBB3',	'ERBB4',	'FGFR1',	'FGFR1',	'FGFR2',	'FGFR2',
                   'FGFR3',	'FGFR3',	'FGFR4',	'FGFR4',	'FLT1',	'FLT4',	'HRAS',	'IGF1R',	'IGF2R',	'IKBKB',
                   'IKBKG',	'INSR',	'KDR',	'KRAS',	'MAP2K1',	'MAP2K2',	'NFKB1',	'NFKB2',	'NGFR',	'NOS3',	'NRAS',
                   'NTRK1',	'NTRK2',	'NTRK3',	'PDGFRA',	'PDGFRB',	'PIK3C2A',	'PIK3C2B',	'PIK3C2G',	'PIK3C3',
                   'PIK3CA',	'PIK3CB',	'PIK3CD',	'PIK3CG',	'PIK3R1',	'PIK3R2',	'PIK3R3',	'PIK3R4',	'PIK3R5',
                   'RPS6KA1',	'RPS6KA2',	'RPS6KA3',	'RPS6KA6',	'TGFBR1',	'TGFBR2',	'TGFBR3',	'TP53',	'YWHAB',
                   'YWHAE',	'YWHAG',	'YWHAH',	'YWHAQ',	'YWHAZ',	'YWHAZ'))
pbmc = Seurat::AddModuleScore(pbmc, features = Apoptosis, name='Apoptosis')
p3<-FeaturePlot(pbmc, features = 'Apoptosis1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10") + 
  ggtitle("Apoptosis", subtitle = paste0("correlation with NES score = ",round(as.numeric(cor.test(pbmc$Apoptosis1,pbmc$NES1)$estimate),2)))
df <- pbmc@meta.data
stat.test <- df %>% t_test(Apoptosis1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
df$seurat_clusters <- factor(df$seurat_clusters, level=c("AC","MES","NPC","OPC"))
p4 <- ggboxplot(df, x = "seurat_clusters", y = "Apoptosis1", fill = "seurat_clusters", 
                palette = col.ls) + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)

Autophagy= list(c('AMBRA1',	'ATG10',	'ATG12',	'ATG13',	'ATG14',	'ATG16L1',
                  'ATG3',	'ATG4A',	'ATG4B',	'ATG5',	'ATG7',	'ATG9A',	'BECN1',
                  'BNIP3L',	'GABARAPL2',	'MAP1LC3A',	'PIK3C3',	'PIK3R4',	'RB1CC1',
                  'SH3GLB1',	'TUBB3',	'ULK1',	'ULK2',	'VMP1',	'WIPI1',	'WIPI2',	'ZFYVE1'))
pbmc = Seurat::AddModuleScore(pbmc, features = Autophagy, name='Autophagy')
p5<-FeaturePlot(pbmc, features = 'Autophagy1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10") + 
  ggtitle("Autophagy", subtitle = paste0("correlation with NES score = ",round(as.numeric(cor.test(pbmc$Autophagy1,pbmc$NES1)$estimate),2)))
df <- pbmc@meta.data
stat.test <- df %>% t_test(Autophagy1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
df$seurat_clusters <- factor(df$seurat_clusters, level=c("AC","MES","NPC","OPC"))
p6 <- ggboxplot(df, x = "seurat_clusters", y = "Autophagy1", fill = "seurat_clusters", 
                palette = col.ls) + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)


Ferroptosis= list(c('ACSL1',	'ACSL3',	'ACSL4',	'ACSL5',	'ACSL6',	'ALOX15',
                    'ATG5',	'ATG7',	'CP',	'CYBB',	'FTH1',	'FTL',	'FTMT',	'GCLC',
                    'GCLM',	'GPX4',	'GSS',	'HMOX1',	'LPCAT3',	'MAP1LC3A',	'MAP1LC3B',
                    'MAP1LC3C',	'NCOA4',	'PCBP1',	'PCBP2',	'PRNP',	'SAT1',	'SAT2',	'SLC11A2',
                    'SLC39A14',	'SLC39A8',	'SLC3A2',	'SLC40A1',	'SLC7A11',	'STEAP3',	'TF',	'TFRC',
                    'TP53',	'VDAC2',	'VDAC3'))
pbmc = Seurat::AddModuleScore(pbmc, features = Ferroptosis, name='Ferroptosis')
p7<-FeaturePlot(pbmc, features = 'Ferroptosis1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")+ 
  ggtitle("Ferroptosis", subtitle = paste0("correlation with NES score = ",round(as.numeric(cor.test(pbmc$Ferroptosis1,pbmc$NES1)$estimate),2)))
df <- pbmc@meta.data
stat.test <- df %>% t_test(Ferroptosis1 ~ seurat_clusters)
stat.test <- stat.test %>% add_xy_position(x = "seurat_clusters")
df$seurat_clusters <- factor(df$seurat_clusters, level=c("AC","MES","NPC","OPC"))
p8 <- ggboxplot(df, x = "seurat_clusters", y = "Ferroptosis1", fill = "seurat_clusters", 
                palette = col.ls) + stat_pvalue_manual (stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.shorten = 0.05)

Necroptosis= list(c('AIFM1',	'ALOX15',	'BAX',	'BCL2',	'BID',	'BIRC2',	'BIRC3',
                    'CAMK2A',	'CAMK2B',	'CAMK2D',	'CAMK2G',	'CAPN1',	'CAPN2',	'CASP1',
                    'CASP8',	'CFLAR',	'CHMP1A',	'CHMP1B',	'CHMP2A',	'CHMP2B',	'CHMP3',	'CHMP4A',
                    'CHMP4B',	'CHMP4C',	'CHMP5',	'CHMP6',	'CHMP7',	'CYBB',	'CYLD',	'DNM1L',	'EIF2AK2',
                    'FADD',	'FAF1',	'FAS',	'FASLG',	'FTH1',	'FTL',	'GLUD1',	'GLUD2',	'GLUL',
                    'H2AB1',	'H2AB2',	'H2AB3',	'H2AC1',	'H2AC11',	'H2AC12',	'H2AC13',	'H2AC14',
                    'H2AC15',	'H2AC16',	'H2AC17',	'H2AC18',	'H2AC19',	'H2AC20',	'H2AC21',	'H2AC4',
                    'H2AC6',	'H2AC7',	'H2AC8',	'H2AJ',	'H2AW',	'H2AX',	'H2AZ1',	'H2AZ2',	'HMGB1',
                    'HSP90AA1',	'HSP90AB1',	'IFNA1',	'IFNA10',	'IFNA13',	'IFNA14',	'IFNA16',	'IFNA17',
                    'IFNA2',	'IFNA21',	'IFNA4',	'IFNA5',	'IFNA6',	'IFNA7',	'IFNA8',	'IFNAR1',
                    'IFNAR2',	'IFNB1',	'IFNG',	'IFNGR1',	'IFNGR2',	'IL1A',	'IL1B',	'IL33',	'IRF9',
                    'JAK1',	'JAK2',	'JAK3',	'JMJD7-PLA2G4B',	'MACROH2A1',	'MACROH2A2',	'MAPK10',
                    'MAPK8',	'MAPK9',	'MLKL',	'NLRP3',	'PARP1',	'PGAM5',	'PLA2G4A',	'PLA2G4B',
                    'PLA2G4C',	'PLA2G4D',	'PLA2G4E',	'PLA2G4F',	'PPIA',	'PPID',	'PYCARD',	'PYGB',
                    'PYGL',	'PYGM',	'RBCK1',	'RIPK1',	'RIPK3',	'RNF103-CHMP3',	'RNF31',	'SHARPIN',
                    'SLC25A31',	'SLC25A4',	'SLC25A5',	'SLC25A6',	'SMPD1',	'SPATA2',	'SPATA2L',
                    'SQSTM1',	'STAT1',	'STAT2',	'STAT3',	'STAT4',	'STAT5A',	'STAT5B',	'STAT6',
                    'TICAM1',	'TICAM2',	'TLR3',	'TLR4',	'TNF',	'TNFAIP3',	'TNFRSF10A',	'TNFRSF10B',
                    'TNFRSF1A',	'TNFSF10',	'TRADD',	'TRAF2',	'TRAF5',	'TRPM7',	'TYK2',	'USP21',
                    'VDAC1',	'VDAC2',	'VDAC3',	'VPS4A',	'VPS4B',	'XIAP',	'ZBP1'))
pbmc = Seurat::AddModuleScore(pbmc, features = Necroptosis, name='Necroptosis')
p9<-FeaturePlot(pbmc, features = 'Necroptosis1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")+ 
  ggtitle("Necroptosis", subtitle = paste0("correlation with NES score = ",round(as.numeric(cor.test(pbmc$Necroptosis1,pbmc$NES1)$estimate),2)))
p10<-VlnPlot(pbmc,features = 'Necroptosis1', group.by = 'seurat_clusters', cols = col.ls, pt.size = 0.05)+ 
  ggtitle("Necroptosis")
cor.test(pbmc$Necroptosis1,pbmc$NES1)

Pyroptosis= list(c('DHX9',	'GSDME',	'GSDMA',	'NLRP9',	'APIP',	'GSDMB',	'GSDMC',	'NLRC4',	'GSDMD',	'AIM2'))
pbmc = Seurat::AddModuleScore(pbmc, features = Pyroptosis, name='Pyroptosis')
p11<-FeaturePlot(pbmc, features = 'Pyroptosis1', reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")+ 
  ggtitle("Pyroptosis", subtitle = paste0("correlation with NES score = ",round(as.numeric(cor.test(pbmc$Pyroptosis1,pbmc$NES1)$estimate),2)))
p12<-VlnPlot(pbmc,features = 'Pyroptosis1', group.by = 'seurat_clusters', cols = col.ls, pt.size = 0.05)+ 
  ggtitle("Pyroptosis")
cor.test(pbmc$Pyroptosis1,pbmc$NES1)

pdf("plot_final/malignant_NES_ferroptosis_umap.pdf", height = 4, width = 15)
FeaturePlot(pbmc, features = c("NES1", "Ferroptosis1"), blend = TRUE, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

saveRDS(pbmc, file = "combined_malignant.rds")

### for BMDM
pdf("plot_new/BMDM_death.pdf")
VlnPlot(pbmc, features = c("Pyroptosis1", "Apoptosis1", "Necroptosis1", "Autophagy1","Ferroptosis1"),
        cols = cluster_color, pt.size = 0)
dev.off()

library(GGally)
data <- data.frame(NES=pbmc$NES1,Apoptosis=pbmc$Apoptosis1,Autophagy=pbmc$Autophagy1,
                   Ferroptosis=pbmc$Ferroptosis1,Necroptosis=pbmc$Necroptosis1,Pyroptosis=pbmc$Pyroptosis1)

library(reshape2)
data2 <- data.frame(NES=pbmc$NES1,Apoptosis=pbmc$Apoptosis1,Autophagy=pbmc$Autophagy1,
                    Ferroptosis=pbmc$Ferroptosis1,Necroptosis=pbmc$Necroptosis1,Pyroptosis=pbmc$Pyroptosis1,
                    group=pbmc$seurat_clusters)
data2<-melt(data2,value.name = "group")
colnames(data2)<-c("cell_type","death_type","value")

pdf("plot_final/malignant_NES_correlation.pdf")
p1+p2
p3+p4
p5+p6
p7+p8
p9+p10
p11+p12
ggcorr(data, method = c("everything", "pearson"), label = TRUE, label_size = 4, label_round = 2, label_alpha = TRUE)
# lapply(unique(pbmc$seurat_cluster),function(x){
#   pos<-which(pbmc$seurat_cluster==x)
#   data <- data.frame(NES=pbmc$NES1[pos],Apoptosis=pbmc$Apoptosis1[pos],Autophagy=pbmc$Autophagy1[pos],
#                      Ferroptosis=pbmc$Ferroptosis1[pos],Necroptosis=pbmc$Necroptosis1[pos],Pyroptosis=pbmc$Pyroptosis1[pos])
#   ggcorr(data, method = c("everything", "pearson"), label = TRUE, label_size = 4, label_round = 2, label_alpha = TRUE)+
#     ggtitle(paste0("seurat cluster ",x))
# })
dev.off()

ggplot(data2, aes(death_type,value,fill=cell_type))+geom_boxplot(alpha=0.5)

M1_M2<-list(M1=c("CCL2","CCL3","CCL4","CCL5","CCL8","CCR7","CD74","CSF2","CXCL10","HLA-DRA","HLA-DRB","IFNG",
                 "IL1B","IL1R1","IL6","INOS","IRF5","NFKB1","TLR2","TLR4","TNF"),
      M2=c("ARG1","CD74","CCL1","CCL17","CCL22","CXCL16","CXCR4","HLA-DRA","HLA-DRB","IL10",
           "IL4","IRF4","MRC1","NFKB1","TGFB1","TNF"))
pbmc = Seurat::AddModuleScore(pbmc, features = M1_M2, name='Macrophage')
FeaturePlot(pbmc, features = c("Macrophage1","Macrophage2"), reduction = 'umap', max.cutoff = "q90", min.cutoff = "q10")
VlnPlot(pbmc, features = c("Macrophage1","Macrophage2"),cols = cluster_color)

pbmc$mac_state <- ifelse(pbmc$Macrophage1>pbmc$Macrophage2,"M1","M2")
pbmc$mac_state[pbmc$Macrophage1<0&pbmc$Macrophage2<0] <- "M0"
DimPlot(pbmc, group.by ="mac_state")

CellInfo <- pbmc@meta.data
P1=CellInfo %>% ggplot(aes(x=sample, fill=mac_state)) +
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P3=CellInfo %>% ggplot(aes(x=grade, fill=mac_state)) +
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P3g=ggplotGrob(P3)
P2=ggplot(CellInfo, aes(mac_state , fill=mac_state))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)
pdf("plot_new/proportion_mac_state.pdf",width=7,height=5)
grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
grid.arrange(grobs=list(P3g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

saveRDS(pbmc, file = "result/BMDM.rds")


