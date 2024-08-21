library(ggplot2)
library(rtracklayer)

# Retrieve chromosome lengths for hg38
hg38_chrom_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
hg38_chrom_lengths <- as.data.frame(hg38_chrom_lengths)
hg38_chrom_lengths$chr <- rownames(hg38_chrom_lengths)

# Filter for main chromosomes (1-22, X, Y)
hg38_chrom_lengths <- hg38_chrom_lengths[grep("^(chr[0-9XY]+)$", rownames(hg38_chrom_lengths)), ]

# Create a cumulative length column for plotting vertical lines
hg38_chrom_lengths$hg38_chrom_lengths <- as.numeric(hg38_chrom_lengths$hg38_chrom_lengths)
hg38_chrom_lengths$cum_length <- cumsum(hg38_chrom_lengths$hg38_chrom_lengths)
hg38_chrom_lengths$cum_start <- c(0, head(hg38_chrom_lengths$cum_length, -1))  # Start position of each chromosome
hg38_chrom_lengths$midpoint <- (hg38_chrom_lengths$cum_start + hg38_chrom_lengths$cum_length) / 2  # Midpoint for labels


cnv_data <- read.csv("data/N20240295-01-01_report/3_Variant/3_CNV/CNV.anno.csv")
cnv_data <- cnv_data[,1:9]

# Adjust CNV data by adding cumulative chromosome lengths
cnv_data$adjusted_start <- cnv_data$Start + hg38_chrom_lengths$cum_start[match(cnv_data$X.Chr, hg38_chrom_lengths$chr)]
cnv_data$adjusted_end <- cnv_data$End + hg38_chrom_lengths$cum_start[match(cnv_data$X.Chr, hg38_chrom_lengths$chr)]

# Create a ggplot object
ggplot(cnv_data, aes(x = adjusted_start, xend = adjusted_end, y = 1, yend = 1, color = CNVType)) + 
  geom_segment(size = 20) +  # Thick lines for CNV regions
  geom_vline(data = hg38_chrom_lengths, aes(xintercept = cum_length), linetype = "dashed", color = "black") +  # Vertical lines
  scale_x_continuous(breaks = hg38_chrom_lengths$midpoint, labels = hg38_chrom_lengths$chr) +  # Chromosome labels
  scale_color_manual(values = c("DUP" = "red", "DEL" = "blue")) +  # Color for DUP
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "Chromosomal Position", y = NULL, color = "CNV Type")+
  facet_wrap(~ SampleName, ncol = 1, strip.position = "left")

