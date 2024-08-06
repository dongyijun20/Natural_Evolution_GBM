library(ROGUE)

nonmalignant <- readRDS("result/nonmalignant_new.rds")
myeloid <- nonmalignant[,nonmalignant$celltype=="Myeloid cells"]

# myeloid <- subset(myeloid, sample%in%c("NJ01_1","NJ01_2","NJ02_1","NJ02_2"))
expr = as.data.frame(myeloid@assays$RNA@counts)
meta <- myeloid@meta.data
ent.res <- SE_fun(expr)
SEplot(ent.res)
table(myeloid$sample)

ggplot(meta, aes(sample , fill=sample))+geom_bar(stat="count",colour = "black",width = 0.7)+  
  scale_fill_manual(values = SampleColor)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"))+labs(y ="number of cells", x= NULL)+ 
  theme(legend.position = "none")


# rogue.value <- CalculateRogue(ent.res, platform = "UMI")
# rogue.value

# sample
rogue.res <- rogue(expr, labels = meta$grade, samples = meta$sample, platform = "UMI", span = 0.6)
rogue.res
rogue.boxplot(rogue.res)

res1 <- rogue.res
res1$sample <- rownames(rogue.res)
res1 <- melt(res1)
res1 <- na.omit(res1)
res1$value <- as.numeric(res1$value)

# subtype
rogue.res <- rogue(expr, labels = meta$subtype, samples = meta$sample, platform = "UMI", span = 0.6)
rogue.res
rogue.boxplot(rogue.res)

res2 <- rogue.res
res2$sample <- rownames(rogue.res)
res2 <- melt(res2)
res2 <- na.omit(res2)
res2$value <- as.numeric(res2$value)

GradeColor <- readRDS("color/GradeColor.rds")
pdf("plot_new/myeloid_rogue_boxplot.pdf")
ggboxplot(res1, x = "variable", y = "value", 
          color = "variable", palette = GradeColor, add = "point")+ ylab("ROGUE")+xlab("Lesions")
ggboxplot(res2, x = "variable", y = "value",
          color = "variable", palette = subtype_color, add = "point")+ ylab("ROGUE")+xlab("Subtypes")
dev.off()
