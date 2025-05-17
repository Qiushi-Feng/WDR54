### Genomic mutation landscape
## CNV data processing
rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(SummarizedExperiment)
library(TCGAbiolinks)

output_dir <- "C:/R/C_Project/F4/CNV"

query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment"
)
GDCdownload(query, method = "api", directory = output_dir)
HNSC_CNV_download <- GDCprepare(
  query = query,
  save = TRUE,
  save.filename = file.path(output_dir, "HNSC_CNV_download.rda"),
  directory = output_dir
)

load(file.path(output_dir, "HNSC_CNV_download.rda"))

tumorCNV <- as.data.frame(HNSC_CNV_download)
tumorCNV <- tumorCNV[, c('Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')]

write.table(
  tumorCNV,
  file = file.path(output_dir, "HNSC_CNV.txt"),
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

group_file_path <- "C:/R/C_Project/F4/CNV/TCGA_GROUP.csv"
cnv_file_path   <- "C:/R/C_Project/F4/CNV/HNSC_CNV.txt"

group_data  <- read.csv(group_file_path)
low_group   <- group_data[, 1]
high_group  <- group_data[, 2]

cnv_data       <- read.table(
  cnv_file_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
cnv_full_names <- cnv_data[, 1]
cnv_short_names <- substr(cnv_full_names, 1, 15)

low_matched <- cnv_data[cnv_short_names %in% low_group, ]
write.table(
  low_matched,
  file = "C:/R/C_Project/F4/CNV/CNV_Low.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

high_matched <- cnv_data[cnv_short_names %in% high_group, ]
write.table(
  high_matched,
  file = "C:/R/C_Project/F4/CNV/CNV_High.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

no_matched <- cnv_data[!cnv_short_names %in% c(low_group, high_group), ]
write.table(
  no_matched,
  file = "C:/R/C_Project/F4/CNV/CNV_no_matched.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("CNV数据匹配和保存完成。")




#随后进行在线GISTIC处理

filename <- "C:/R/C_Project/F4/CNV/GISTIC/snp6.na35.remap.hg38.subset.txt"
finalResultName <- "C:/R/C_Project/F4/CNV/GISTIC/marker_file.txt"

read_file <- file(filename, open = "r")
out_file <- file(finalResultName, open = "w")

while (length(line <- readLines(read_file, n = 1, warn = FALSE)) > 0) {
  data <- strsplit(line, "\\s+")[[1]]
  if (data[6] == "FALSE") {
    writeLines(paste(data[1], data[2], data[3], sep = "\t"), out_file)
  }
}

close(read_file)
close(out_file)

file_path <- "C:/R/C_Project/F4/CNV/GISTIC/marker_file.txt"
marker_data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(marker_data) <- c("Maker Name", "Chromosome", "Maker Position")
write.table(marker_data, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



#随后进行在线GENEPATTERN处理

rm(list = ls())
options(stringsAsFactors = F)
library(maftools)

gistic_dir <- "C:/R/C_Project/F4/CNV/GISTIC/CNV_Low/600219"
laml.gistic <- readGistic(
  gisticAllLesionsFile = file.path(gistic_dir, "all_lesions.conf_90.txt"),
  gisticAmpGenesFile   = file.path(gistic_dir, "amp_genes.conf_90.txt"),
  gisticDelGenesFile   = file.path(gistic_dir, "del_genes.conf_90.txt"),
  gisticScoresFile     = file.path(gistic_dir, "scores.gistic"),
  isTCGA               = TRUE
)
output_dir <- gistic_dir

original_width   <- 800
original_height  <- 800
new_width_other  <- original_width * 2/3
new_height_other <- original_height * 1/3
new_width_bubble <- original_width * 2/3
new_height_bubble<- original_height * 1/2

png(file.path(output_dir, "gisticChromPlot.png"), width = new_width_other, height = new_height_other)
gisticChromPlot(gistic = laml.gistic, ref.build = "hg38")
dev.off()

png(file.path(output_dir, "gisticBubblePlot.png"), width = new_width_bubble, height = new_height_bubble)
gisticBubblePlot(gistic = laml.gistic)
dev.off()

png(file.path(output_dir, "gisticOncoPlot.png"), width = new_width_other, height = new_height_other)
gisticOncoPlot(gistic = laml.gistic, sortByAnnotation = TRUE, top = 10)
dev.off()

segment_file_path <- file.path(gistic_dir, "CNV_Low.txt")
if (file.exists(segment_file_path)) {
  png(file.path(output_dir, "plotCBSsegments.png"), width = new_width_other, height = new_height_other)
  plotCBSsegments(cbsFile = segment_file_path, tsb = "ALL", ref.build = "hg38")
  dev.off()
} else {
  cat("Error: CNV_Low.txt 文件不存在或不可读。")
}

rm(list = ls())
options(stringsAsFactors = F)
library(maftools)

gistic_dir <- "C:/R/C_Project/F4/CNV/GISTIC/CNV-High/600218"
laml.gistic <- readGistic(
  gisticAllLesionsFile = file.path(gistic_dir, "all_lesions.conf_90.txt"),
  gisticAmpGenesFile   = file.path(gistic_dir, "amp_genes.conf_90.txt"),
  gisticDelGenesFile   = file.path(gistic_dir, "del_genes.conf_90.txt"),
  gisticScoresFile     = file.path(gistic_dir, "scores.gistic"),
  isTCGA               = TRUE
)
output_dir <- gistic_dir

original_width   <- 800
original_height  <- 800
new_width_other  <- original_width * 2/3
new_height_other <- original_height * 1/3
new_width_bubble <- original_width * 2/3
new_height_bubble<- original_height * 1/2

png(file.path(output_dir, "gisticChromPlot.png"), width = new_width_other, height = new_height_other)
gisticChromPlot(gistic = laml.gistic, ref.build = "hg38")
dev.off()

png(file.path(output_dir, "gisticBubblePlot.png"), width = new_width_bubble, height = new_height_bubble)
gisticBubblePlot(gistic = laml.gistic)
dev.off()

png(file.path(output_dir, "gisticOncoPlot.png"), width = new_width_other, height = new_height_other)
gisticOncoPlot(gistic = laml.gistic, sortByAnnotation = TRUE, top = 10)
dev.off()

segment_file_path <- file.path(gistic_dir, "CNV_High.txt")
if (file.exists(segment_file_path)) {
  png(file.path(output_dir, "plotCBSsegments.png"), width = new_width_other, height = new_height_other)
  plotCBSsegments(cbsFile = segment_file_path, tsb = "ALL", ref.build = "hg38")
  dev.off()
} else {
  cat("Error: CNV_High.txt 文件不存在或不可读。")
}

## Mutation visualization
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
library(ggplot2)

file_path <- "C:/R/C_Project/F4/UALCAN/promoter_methylation_Normal_VS_Tumor.csv"
data <- read.csv(file_path)
names(data)
colors <- c(
  "Primary tumor<br>(n=528)" = rgb(222, 82, 108, maxColorValue = 255),
  "Normal<br>(n=50)"           = rgb(126, 166, 217, maxColorValue = 255)
)
output_file   <- "C:/R/C_Project/F4/UALCAN/promoter_methylation_boxplot.png"
current_width <- 8
new_width     <- current_width * 1/2

p <- ggplot(
  data,
  aes(
    x      = factor(TCGA.samples, levels = data$TCGA.samples),
    ymin   = Series.1..low.,
    lower  = Series.1..q1.,
    middle = Series.1..median.,
    upper  = Series.1..q3.,
    ymax   = Series.1..high.,
    fill   = TCGA.samples
  )
) +
  geom_boxplot(stat = "identity", color = "black") +
  labs(x = "Sample Type", y = "Beta Value") +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_blank(),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position  = "none"
  )

ggsave(filename = output_file, plot = p, width = new_width, height = 6)
cat("图形已保存为: ", output_file)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
library(ggplot2)

file_path <- "C:/R/C_Project/F4/UALCAN/WDR54_HPV.csv"
data <- read.csv(file_path)
names(data)
colors <- c(
  "HPV+ve<br>(n=80)" = rgb(248, 242, 164, maxColorValue = 255),
  "HPV-ve<br>(n=434)" = rgb(212, 230, 188, maxColorValue = 255),
  "Normal<br>(n=44)"  = rgb(126, 166, 217, maxColorValue = 255)
)
output_file   <- "C:/R/C_Project/F4/UALCAN/WDR54_HPV_boxplot.png"
current_width <- 10
new_width     <- current_width * 1/2

p <- ggplot(
  data,
  aes(
    x      = factor(`TCGA.samples`, levels = data$`TCGA.samples`),
    ymin   = `Series.1..low.`,
    lower  = `Series.1..q1.`,
    middle = `Series.1..median.`,
    upper  = `Series.1..q3.`,
    ymax   = `Series.1..high.`,
    fill   = `TCGA.samples`
  )
) +
  geom_boxplot(stat = "identity", color = "black") +
  labs(x = "Sample Type", y = "Transcript per Million") +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_blank(),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position  = "none"
  )

ggsave(filename = output_file, plot = p, width = new_width, height = 6)
cat("图形已保存为: ", output_file)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
library(ggplot2)

file_path <- "C:/R/C_Project/F4/UALCAN/WDR54_TP53.csv"
data <- read.csv(file_path)
names(data)
colors <- c(
  "Normal<br>(n=44)"             = rgb(126, 166, 217, maxColorValue = 255),
  "TP53-Mutant<br>(n=327)"        = rgb(181, 140, 154, maxColorValue = 255),
  "TP53-NonMutant<br>(n=175)"     = rgb(137, 133, 183, maxColorValue = 255)
)
output_file   <- "C:/R/C_Project/F4/UALCAN/WDR54_TP53_boxplot.png"
current_width <- 10
new_width     <- current_width * 1/2

p <- ggplot(
  data,
  aes(
    x      = factor(`TCGA.samples`, levels = data$`TCGA.samples`),
    ymin   = `Series.1..low.`,
    lower  = `Series.1..q1.`,
    middle = `Series.1..median.`,
    upper  = `Series.1..q3.`,
    ymax   = `Series.1..high.`,
    fill   = `TCGA.samples`
  )
) +
  geom_boxplot(stat = "identity", color = "black") +
  labs(x = "Sample Type", y = "Transcript per Million") +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_blank(),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position  = "none"
  )

ggsave(filename = output_file, plot = p, width = new_width, height = 6)
cat("图形已保存为: ", output_file)

##SNV analysis
library(TCGAbiolinks)

query <- GDCquery(
    project = "TCGA-HNSC", 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)
GDCdownload(query, directory = "C:/R/C_Project/F4/TCGA_genetic_mutation")
maf_data <- GDCprepare(
    query, 
    save = TRUE, 
    save.filename = "C:/R/C_Project/F4/TCGA_genetic_mutation/TCGA-HNSC_SNP.Rdata",
    directory = "C:/R/C_Project/F4/TCGA_genetic_mutation"
)
csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/TCGA-HNSC_SNP.csv"
write.csv(maf_data, file = csv_file_path, row.names = FALSE)
cat("数据已保存为 RData 和 CSV 格式。\n")

library(openxlsx)
library(tidyverse)
library(limma)
library(readr)

data_dir <- "C:/R/C_Project/F4/TCGA_genetic_mutation"
TCGA_rawdata <- read_tsv(file.path(data_dir, "TCGA-HNSC.htseq_counts.tsv"))
dim(TCGA_rawdata)
probeMap <- read.table(
    file.path(data_dir, "gencode.v22.annotation.gene.probeMap"),
    sep = "\t", header = TRUE
)
TCGA_gset <- TCGA_rawdata %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene, starts_with("TCGA"))
TCGA_gset <- as.data.frame(avereps(TCGA_gset[,-1], ID = TCGA_gset$gene))
colnames(TCGA_gset) <- substring(colnames(TCGA_gset), 1, 15) %>% gsub("-", ".", .)
output_file_path <- file.path(data_dir, "TCGA_HNSC_Countdata_log2+1.csv")
write.csv(TCGA_gset, file = output_file_path, row.names = TRUE)
TCGA_gset[1:4, 1:4]
TCGA_group_list <- ifelse(
  as.numeric(substring(colnames(TCGA_gset), 14, 15)) < 10,
  "Tumor", "Normal"
) %>% factor(., levels = c("Normal", "Tumor"))
table(TCGA_group_list)
cat("处理后的数据已保存为:", output_file_path, "\n")

csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/TCGA_GROUP.csv"
data <- read.csv(csv_file_path)
sample_names <- data[[1]]
sample_types <- sapply(strsplit(sample_names, "\\."), function(x) x[4])
filtered_data <- data[sample_types == "01", ]
write.csv(filtered_data, file = csv_file_path, row.names = FALSE)
cat("已删除正常组织样本，结果已覆盖原文件:", csv_file_path, "\n")

group_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/TCGA_GROUP.csv"
snp_file_path   <- "C:/R/C_Project/F4/TCGA_genetic_mutation/TCGA-HNSC_SNP.csv"
group_data  <- read.csv(group_file_path)
low_group   <- group_data[, 1]
high_group  <- group_data[, 2]
snp_data    <- read.csv(snp_file_path)
snp_full_names <- snp_data[, 16]
snp_short_names <- substr(snp_full_names, 1, 15)
low_matched  <- snp_data[snp_short_names %in% low_group, ]
write.csv(low_matched, file = "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv", row.names = FALSE)
high_matched <- snp_data[snp_short_names %in% high_group, ]
write.csv(high_matched, file = "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv", row.names = FALSE)
no_matched   <- snp_data[!snp_short_names %in% c(low_group, high_group), ]
write.csv(no_matched, file = "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_no_matched.csv", row.names = FALSE)
cat("匹配和保存完成。\n")

library(maftools)

csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
high_maf    <- read.csv(csv_file_path)
maf_object  <- maftools::read.maf(maf = high_maf)
png_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High_oncoplot.png"
png(filename = png_file_path, width = 1000, height = 800)
maftools::oncoplot(maf = maf_object)
dev.off()
cat("Oncoplot图形已保存为：", png_file_path, "\n")

csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
low_maf     <- read.csv(csv_file_path)
maf_object  <- maftools::read.maf(maf = low_maf)
png_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low_oncoplot.png"
png(filename = png_file_path, width = 1000, height = 800)
maftools::oncoplot(maf = maf_object)
dev.off()
cat("Oncoplot图形已保存为：", png_file_path, "\n")

library(maftools)
high_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
low_csv_file_path  <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_high       <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_High_summary.png"
output_low        <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_summary.png"
high_maf          <- read.csv(high_csv_file_path)
high_maf_object   <- maftools::read.maf(maf = high_maf)
png(filename = output_high, width = 1200, height = 800)
plotmafSummary(maf = high_maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
cat("SNV_High的突变数据总结图已保存为：", output_high, "\n")
low_maf         <- read.csv(low_csv_file_path)
low_maf_object  <- maftools::read.maf(maf = low_maf)
png(filename = output_low, width = 1200, height = 800)
plotmafSummary(maf = low_maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
cat("SNV_Low的突变数据总结图已保存为：", output_low, "\n")

library(maftools)
high_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
low_csv_file_path  <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_high_tiTv   <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_High_TiTv.png"
output_low_tiTv    <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_TiTv.png"
high_maf           <- read.csv(high_csv_file_path)
high_maf_object    <- maftools::read.maf(maf = high_maf)
high_titv          <- titv(maf = high_maf_object, plot = FALSE, useSyn = TRUE)
png(filename = output_high_tiTv, width = 1200, height = 800)
plotTiTv(res = high_titv)
dev.off()
cat("SNV_High的Ti/Tv比率图已保存为：", output_high_tiTv, "\n")
low_maf            <- read.csv(low_csv_file_path)
low_maf_object     <- maftools::read.maf(maf = low_maf)
low_titv           <- titv(maf = low_maf_object, plot = FALSE, useSyn = TRUE)
png(filename = output_low_tiTv, width = 1200, height = 800)
plotTiTv(res = low_titv)
dev.off()
cat("SNV_Low的Ti/Tv比率图已保存为：", output_low_tiTv, "\n")
library(maftools)
high_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
low_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_file_path_high <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_High_rainfall.png"
output_file_path_low <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_rainfall.png"

high_maf <- read.csv(high_csv_file_path)
high_maf_object <- maftools::read.maf(maf = high_maf)
png(filename = output_file_path_high, width = 1200, height = 800)
rainfallPlot(maf = high_maf_object, detectChangePoints = TRUE, pointSize = 0.6)
dev.off()
cat("SNV_High的雨量图已保存为：", output_file_path_high, "\n")

low_maf <- read.csv(low_csv_file_path)
low_maf_object <- maftools::read.maf(maf = low_maf)
png(filename = output_file_path_low, width = 1200, height = 800)
rainfallPlot(maf = low_maf_object, detectChangePoints = TRUE, pointSize = 0.6)
dev.off()
cat("SNV_Low的雨量图已保存为：", output_file_path_low, "\n")

if (!requireNamespace("wordcloud", quietly = TRUE)) {
    install.packages("wordcloud")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
}
library(wordcloud)
library(RColorBrewer)
library(maftools)
library(wordcloud)

high_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
low_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_file_path_high <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_High_geneCloud.png"
output_file_path_low <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_geneCloud.png"

high_maf <- read.csv(high_csv_file_path)
low_maf <- read.csv(low_csv_file_path)
high_maf_object <- maftools::read.maf(maf = high_maf)
low_maf_object <- maftools::read.maf(maf = low_maf)

high_genes <- getGeneSummary(high_maf_object)
low_genes <- getGeneSummary(low_maf_object)
high_genes_filtered <- high_genes[high_genes$total >= 15, ]
low_genes_filtered <- low_genes[low_genes$total >= 15, ]

png(filename = output_file_path_high, width = 1200, height = 800)
wordcloud(words = high_genes_filtered$Hugo_Symbol, freq = high_genes_filtered$total, min.freq = 15, colors = brewer.pal(8, "Dark2"))
dev.off()

png(filename = output_file_path_low, width = 1200, height = 800)
wordcloud(words = low_genes_filtered$Hugo_Symbol, freq = low_genes_filtered$total, min.freq = 15, colors = brewer.pal(8, "Dark2"))
dev.off()

cat("HNSC_WDR54_High和HNSC_WDR54_Low的基因云图已分别保存为：\n", output_file_path_high, "\n", output_file_path_low, "\n")

Our_maf <- read.csv("Our_maf.csv", header = TRUE)
our_maf <- read.maf(maf = Our_maf)
pt.vs.rt <- mafCompare(m1 = laml, m2 = our_maf, m1Name = 'LIHC', m2Name = 'OUR', minMut = 5)
print(pt.vs.rt)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.01, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

library(maftools)
csv_file_path_primary <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
csv_file_path_relapse <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_dir <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary"

primary_maf <- read.csv(csv_file_path_primary, stringsAsFactors = FALSE)
relapse_maf <- read.csv(csv_file_path_relapse, stringsAsFactors = FALSE)
primary_maf_object <- maftools::read.maf(maf = primary_maf)
relapse_maf_object <- maftools::read.maf(maf = relapse_maf)

pt_vs_rt <- mafCompare(m1 = primary_maf_object, m2 = relapse_maf_object, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt_vs_rt)

output_file <- file.path(output_dir, "Primary_vs_Relapse_CoBarplot.png")
png(filename = output_file, width = 1200, height = 800)
coBarplot(m1 = primary_maf_object, m2 = relapse_maf_object, m1Name = "Primary", m2Name = "Relapse")
dev.off()
cat("条形图已成功保存到:", output_file)

library(maftools)
library(data.table)

high_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
low_csv_file_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
output_csv_path_high <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_High_OncogenicPathways_matrix.csv"
output_csv_path_low <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_OncogenicPathways_matrix.csv"

high_maf <- read.csv(high_csv_file_path)
low_maf <- read.csv(low_csv_file_path)
high_maf_object <- read.maf(maf = high_maf)
low_maf_object <- read.maf(maf = low_maf)

high_pathways <- pathways(maf = high_maf_object)
low_pathways <- pathways(maf = low_maf_object)
high_pathways_matrix <- as.matrix(high_pathways)
low_pathways_matrix <- as.matrix(low_pathways)

write.csv(high_pathways_matrix, file = output_csv_path_high, row.names = FALSE)
write.csv(low_pathways_matrix, file = output_csv_path_low, row.names = FALSE)

cat("SNV_High的致癌信号通路分析结果已保存为矩阵格式的CSV文件：", output_csv_path_high, "\n")
cat("SNV_Low的致癌信号通路分析结果已保存为矩阵格式的CSV文件：", output_csv_path_low, "\n")
library(ggplot2)

input_csv_path <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_OncogenicPathways_matrix.csv"
output_png_path_with_labels <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/Affected_Genes_vs_Pathway_Size_with_labels.png"
output_png_path_without_labels <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/Affected_Genes_vs_Pathway_Size_without_labels.png"

data <- read.csv(input_csv_path)
num_points <- nrow(data[!is.na(data$Pathway), ])
cat("散点图中总共有", num_points, "个散点。\n")

p_with_labels <- ggplot(data, aes(x = N, y = n_affected_genes, color = Group, size = Fraction_mutated_samples, label = Pathway)) +
  geom_point(alpha = 0.7) +
  geom_text(vjust = -1, hjust = 1, size = 3) +
  labs(title = "Affected_Genes vs Pathway_Size",
       x = "Pathway_Size (N)",
       y = "Affected_Genes (n_affected_genes)") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal()

ggsave(output_png_path_with_labels, plot = p_with_labels, width = 12, height = 7, dpi = 300)

p_without_labels <- ggplot(data, aes(x = N, y = n_affected_genes, color = Group, size = Fraction_mutated_samples)) +
  geom_point(alpha = 0.7) +
  labs(title = "Affected_Genes vs Pathway_Size",
       x = "Pathway_Size (N)",
       y = "Affected_Genes (n_affected_genes)") +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal()

ggsave(output_png_path_without_labels, plot = p_without_labels, width = 12, height = 7, dpi = 300)

cat("带有Pathway名称的散点图已保存为PNG文件：", output_png_path_with_labels, "\n")
cat("不带有Pathway名称的散点图已保存为PNG文件：", output_png_path_without_labels, "\n")

library(maftools)

csv_file_path_low <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_Low.csv"
snp_data_low     <- read.csv(csv_file_path_low)
maf_object_low   <- read.maf(maf = snp_data_low)

csv_file_path_high <- "C:/R/C_Project/F4/TCGA_genetic_mutation/SNV_High.csv"
snp_data_high     <- read.csv(csv_file_path_high)
maf_object_high   <- read.maf(maf = snp_data_high)

gene_of_interest <- "TP53"
output_png_path  <- "C:/R/C_Project/F4/TCGA_genetic_mutation/Supplementary/SNV_Low_High_TP53_lollipop.png"

png(output_png_path, width = 1200, height = 800)
lollipopPlot2(
  m1      = maf_object_low,
  m2      = maf_object_high,
  gene    = gene_of_interest,
  AACol1  = 'HGVSp_Short',
  AACol2  = 'HGVSp_Short',
  m1_name = "SNV_Low",
  m2_name = "SNV_High"
)
dev.off()














