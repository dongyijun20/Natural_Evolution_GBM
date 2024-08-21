
malignant <- readRDS("result/malignant_new.rds")
nonmalignant <- readRDS("result/nonmalignant_new.rds")
myeloid <- readRDS("result/myeloid.rds")

sce.markers <- read.table("result/BMDM_markers200.txt", header = T)
library(clusterProfiler)
library(org.Hs.eg.db)

ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
#View(sce.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG", 
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.y = element_text(size=8),axis.text.x = element_text(
  vjust = 0.5, hjust = 0.5
))

Idents(Mo) <- "seurat_clusters"
expr <- AverageExpression(Mo, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
head(expr)

library(msigdbr)
msigdbr_species() #列出有的物种

#选择基因集合
human_KEGG = msigdbr(species = "Homo sapiens", #物种
                     category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol或者ID
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)#list
length(human_KEGG_Set)
head(human_KEGG_Set)

# human_KEGG <- read.gmt("reference/c2.cp.kegg.v7.5.1.symbols.gmt")
# human_KEGG$term <- as.factor(human_KEGG$term)
# human_KEGG_Set <- lapply(unique(human_KEGG$term), function(x){print(x);human_KEGG$gene[human_KEGG$term == x]})
# names(human_KEGG_Set) <- unique(human_KEGG$term)
# 
# length(human_KEGG_Set)
# head(human_KEGG_Set)

library(GSVA)
gsva.kegg <- gsva(expr, gset.idx.list = human_KEGG_Set, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=1)
head(gsva.kegg)

library(pheatmap)
pdf("plot_new/Mo_GSVA_Hallmark_younger.pdf", height = 5, width = 8)
pheatmap(gsva.kegg, show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

## limma------
library(limma)
# 构建分组文件
group_list <- data.frame(celltype = colnames(gsva.kegg), group = colnames(gsva.kegg))

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva.kegg)

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(older-younger, levels = design)

# 差异分析，older vs. younger
fit <- lmFit(gsva.kegg, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

#设置分组
diff$group <- ifelse( diff$logFC > 0 & diff$P.Value < 0.01 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$P.Value < 0.01 ,"down","noSig")
)

diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>% #注释掉该行即可
  arrange(t) %>% 
  rownames_to_column("ID")

diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x =0,nudge_y =0,hjust =diff2$hjust,
            size = 3)+ #设置字体大小
  labs(x = "HALLMARK pathways", #设置标题 和 坐标
       y=paste0("t value of GSVA score\n","celltype-unknown"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none" #去掉legend
  )


