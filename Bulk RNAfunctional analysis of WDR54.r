###pathwaynetwork
# 第一步——数据准备
data <- read.csv("C:/R/C_Project/F5/Network_Heatmap/CCLE_HNSCC_ESCC.csv", row.names = 1)
transposed_data <- t(data)
write.csv(transposed_data, "C:/R/C_Project/F5/Network_Heatmap/CCLE_HNSCC_ESCC.csv", row.names = TRUE)

# 第二步——通路激活计算（CCLE）
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("progeny")
library(progeny)
gene_expression <- as.matrix(read.csv("C:/R/C_Project/F5/Network_Heatmap/CCLE_HNSCC_ESCC.csv", row.names = 1))
pathways <- progeny(gene_expression, scale = TRUE, organism = "Human", top = 100, perm = 1)
controls <- rep(c(TRUE, FALSE), 4)
ctl_mean <- apply(pathways[controls, ], 2, mean)
ctl_sd <- apply(pathways[controls, ], 2, sd)
pathways <- t(apply(pathways, 1, function(x) x - ctl_mean))
pathways <- apply(pathways, 1, function(x) x / ctl_sd)
output_path <- "C:/R/C_Project/F5/Network_Heatmap/CCLE_Pathways_Results.csv"
write.csv(pathways, file = output_path, row.names = TRUE)
print(paste("计算结果已保存至:", output_path))

# 第二步——通路激活计算（GSE65858）
library(progeny)
input_file <- "C:/R/C_Project/F5/Network_Heatmap/GSE65858_clean.csv"
gene_expression <- as.matrix(read.csv(input_file, row.names = 1))
print("基因表达数据的前几行：")
print(head(gene_expression))
pathways <- progeny(gene_expression, scale = TRUE, organism = "Human", top = 100, perm = 1)
print("通路活性计算结果的前几行和前几列：")
print(head(pathways)[1:5, 1:5])
print(paste("计算得到的通路数量：", nrow(pathways)))
controls <- rep(c(TRUE, FALSE), 4)
ctl_mean <- apply(pathways[controls, ], 2, mean)
ctl_sd <- apply(pathways[controls, ], 2, sd)
print("控制组的均值：")
print(ctl_mean)
print("控制组的标准差：")
print(ctl_sd)
pathways_normalized <- t(apply(pathways, 1, function(x) x - ctl_mean))
pathways_normalized <- apply(pathways_normalized, 1, function(x) x / ctl_sd)
pathways_normalized_df <- as.data.frame(pathways_normalized)
print("标准化后的通路活性数据的前几行和前几列：")
print(head(pathways_normalized_df)[1:5, 1:5])
output_file <- "C:/R/C_Project/F5/Network_Heatmap/GSE65858_Pathways_Results.csv"
write.csv(pathways_normalized_df, file = output_file, row.names = TRUE)
print(paste("计算结果已保存至:", output_file))

# 第二步——通路激活计算（TCGA）
library(progeny)
gene_expression <- as.matrix(read.csv("C:/R/C_Project/F5/Network_Heatmap/TCGA_HNSC_mRNA_clean.csv", row.names = 1))
pathways <- progeny(gene_expression, scale = TRUE, organism = "Human", top = 100, perm = 1)
controls <- rep(c(TRUE, FALSE), 4)
ctl_mean <- apply(pathways[controls, ], 2, mean)
ctl_sd <- apply(pathways[controls, ], 2, sd)
pathways <- t(apply(pathways, 1, function(x) x - ctl_mean))
pathways <- apply(pathways, 1, function(x) x / ctl_sd)
output_path <- "C:/R/C_Project/F5/Network_Heatmap/TCGA_Pathways_Results.csv"
write.csv(pathways, file = output_path, row.names = TRUE)
print(paste("计算结果已保存至:", output_path))

# 第三步——数据处理（仅留下代码框架，按需补充）
# （请在此处插入转置结果并添加 E_SCORE, M_SCORE, EMT_SCORE, WDR54 注释的代码）

# 第四步——可视化处理
install.packages("vegan")
install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
library(vegan)
library(devtools)
library(linkET)
library(ggplot2)
library(dplyr)

process_and_plot <- function(file_path, output_csv_path, output_plot_path) {
  data <- read.csv(file_path, sep = ",", header = TRUE)
  pathways_data <- data[, 2:15]
  features_data <- data[, 16:19]
  pearson_results <- data.frame()
  for (i in 1:ncol(features_data)) {
    for (j in 1:ncol(pathways_data)) {
      cor_test <- cor.test(features_data[, i], pathways_data[, j])
      pearson_results <- rbind(
        pearson_results,
        data.frame(
          Feature = colnames(features_data)[i],
          Pathway = colnames(pathways_data)[j],
          R = cor_test$estimate,
          P = cor_test$p.value
        )
      )
    }
  }
  pearson_results <- pearson_results %>% filter(P < 0.05)
  pearson_results <- pearson_results %>%
    mutate(
      pd = cut(P, breaks = c(0.05, 0.01, 0.001, -Inf), labels = c("0.01 - 0.05", "0.001 - 0.01", "< 0.001")),
      rd = cut(abs(R), breaks = c(-Inf, 0.3, 0.5, Inf), labels = c("Weak", "Moderate", "Strong"))
    )
  write.csv(pearson_results, file = output_csv_path, row.names = FALSE)
  custom_colors <- c("0.01 - 0.05" = "#549F9A", "0.001 - 0.01" = "#FD763F", "< 0.001" = "#C30078")
  heatmap_plot <- qcorrplot(correlate(pathways_data), type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(aes(colour = pd, size = rd), data = pearson_results, curvature = nice_curvature()) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = custom_colors) +
    guides(
      size = guide_legend(title = "Correlation strength", override.aes = list(colour = "grey35"), order = 2),
      colour = guide_legend(title = "p-value", override.aes = list(size = 3), order = 1),
      fill = guide_colorbar(title = "Pearson's r", order = 3)
    )
  ggsave(output_plot_path, plot = heatmap_plot, width = 21, height = 16)
}

process_and_plot(
  file_path = "C:/R/C_Project/F5/Network_Heatmap/CCLE_Pathways_Results.csv",
  output_csv_path = "C:/R/C_Project/F5/Network_Heatmap/CCLE_Pearson_Results.csv",
  output_plot_path = "C:/R/C_Project/F5/Network_Heatmap/CCLE_Correlation_Heatmap.png"
)
process_and_plot(
  file_path = "C:/R/C_Project/F5/Network_Heatmap/GSE65858_Pathways_Results.csv",
  output_csv_path = "C:/R/C_Project/F5/Network_Heatmap/GSE65858_Pearson_Results.csv",
  output_plot_path = "C:/R/C_Project/F5/Network_Heatmap/GSE65858_Correlation_Heatmap.png"
)
process_and_plot(
  file_path = "C:/R/C_Project/F5/Network_Heatmap/TCGA_Pathways_Results.csv",
  output_csv_path = "C:/R/C_Project/F5/Network_Heatmap/TCGA_Pearson_Results.csv",
  output_plot_path = "C:/R/C_Project/F5/Network_Heatmap/TCGA_Correlation_Heatmap.png"
)

## relation heatmap





library(dplyr)
input_file <- "C:/R/C_Project/F5/Correlation_Heatmap/gene_list.csv"
output_dir <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA"
gene_list <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
epithelial_genes <- gene_list %>% filter(`Phenotype..marker.` == "Epithelial") %>% select(Symbol) %>% unlist()
mesenchymal_genes <- gene_list %>% filter(`Phenotype..marker.` == "Mesenchymal") %>% select(Symbol) %>% unlist()
write_gmt <- function(gene_list, gene_set_name, output_file) {
  gmt_line <- paste(c(gene_set_name, "NA", gene_list), collapse = "\t")
  writeLines(gmt_line, output_file)
}
epithelial_output_file <- file.path(output_dir, "Epithelial_genes.gmt")
write_gmt(epithelial_genes, "Epithelial_Gene_Set", epithelial_output_file)
mesenchymal_output_file <- file.path(output_dir, "Mesenchymal_genes.gmt")
write_gmt(mesenchymal_genes, "Mesenchymal_Gene_Set", mesenchymal_output_file)
cat("GMT文件已成功生成并保存到:", output_dir)

library(clusterProfiler)
library(GSVA)
library(data.table)
base_path <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/"
expr_files <- list("TCGA.csv", "CCLE.csv", "GSE65858.csv")
gene_set_files <- list("Mesenchymal_genes.gmt", "Epithelial_genes.gmt", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt")
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
run_ssgsea <- function(expr_file, gene_set_file) {
  gene_set_path <- paste0(base_path, gene_set_file)
  if (grepl("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_set_file)) {
    lines <- readLines(gene_set_path)
    if (length(lines)>0 && !nzchar(tail(lines,1))) writeLines(lines, gene_set_path) else writeLines(c(lines,""), gene_set_path)
  }
  genesets <- clusterProfiler::read.gmt(gene_set_path)
  genesets4gsva <- split(genesets$gene, genesets$term)
  expr <- fread(paste0(base_path, expr_file), header=TRUE, data.table=FALSE, check.names=FALSE)
  rownames(expr) <- expr[,1]
  expr[,-1] <- apply(expr[,-1], 2, function(x) as.numeric(as.character(x)))
  expr_matrix <- as.matrix(expr[,-1])
  ssgsea_param <- ssgseaParam(expr_matrix, genesets4gsva, alpha=0.25, normalize=TRUE)
  ssgsea_res <- gsva(ssgsea_param)
  expr_geneset <- as.data.frame(ssgsea_res)
  output_file <- paste0(output_path, gsub(".csv","",expr_file), "_", gsub(".gmt","",gene_set_file), "_ESCORE.csv")
  write.csv(expr_geneset, output_file, row.names=TRUE)
  return(output_file)
}
results <- list()
for (expr_file in expr_files) {
  for (gene_set_file in gene_set_files) {
    result_file <- run_ssgsea(expr_file, gene_set_file)
    results <- c(results, result_file)
    print(paste("结果已保存到：", result_file))
  }
}
print("所有结果文件：")
print(results)

library(data.table)
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
file_paths <- list(
  "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/CCLE_Result.csv",
  "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/TCGA_Result.csv",
  "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/GSE65858_Result.csv"
)
process_file <- function(file_path) {
  data <- read.csv(file_path, header=TRUE)
  data_normalized <- data
  data_normalized[,2:4] <- lapply(data[,2:4], normalize)
  final_data <- data.frame(data_normalized[,1], data_normalized[,2:4], data_normalized[,5])
  colnames(final_data) <- colnames(data)[1:5]
  write.csv(final_data, file=file_path, row.names=FALSE)
  cat("处理完成，已将归一化数据写回原文件：", file_path, "\n")
}
for (file_path in file_paths) process_file(file_path)

library(ggpubr)
calculate_correlation <- function(x_var, y_var, data) {
  pearson_test <- cor.test(data[[x_var]], data[[y_var]], method="pearson")
  spearman_test <- cor.test(data[[x_var]], data[[y_var]], method="spearman")
  data.frame(
    Variable=y_var,
    Pearson_R=pearson_test$estimate, Pearson_P_value=pearson_test$p.value,
    Spearman_R=spearman_test$estimate, Spearman_P_value=spearman_test$p.value
  )
}
plot_and_save <- function(data, x_var, y_var, color, output_file) {
  ggsave(output_file,
         plot=ggscatter(data, x=x_var, y=y_var, add="reg.line", conf.int=TRUE,
                        color=color, size=3, alpha=0.6, ggtheme=theme_minimal()),
         width=9, height=6)
}
process_dataset <- function(data_file, output_prefix, color_escore, color_mscore, color_emtscore) {
  data <- read.csv(data_file)
  result_e <- calculate_correlation("WDR54", "E_SCORE", data)
  result_m <- calculate_correlation("WDR54", "M_SCORE", data)
  result_emt <- calculate_correlation("WDR54", "EMT_SCORE", data)
  final_result <- rbind(result_e, result_m, result_emt)
  write.csv(final_result, file=paste0(output_prefix, "_Correlation_Results.csv"), row.names=FALSE)
  plot_and_save(data, "WDR54", "E_SCORE", color_escore, paste0(output_prefix, "_ESCORE.png"))
  plot_and_save(data, "WDR54", "M_SCORE", color_mscore, paste0(output_prefix, "_MSCORE.png"))
  plot_and_save(data, "WDR54", "EMT_SCORE", color_emtscore, paste0(output_prefix, "_EMTSCORE.png"))
}
datasets <- list(
  list(file="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/CCLE_Result.csv",
       prefix="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/CCLE",
       escore_color=rgb(148,181,216,255),
       mscore_color=rgb(213,123,112,255),
       emtscore_color=rgb(172,210,199,255)),
  list(file="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/TCGA_Result.csv",
       prefix="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/TCGA",
       escore_color=rgb(148,181,216,255),
       mscore_color=rgb(213,123,112,255),
       emtscore_color=rgb(172,210,199,255)),
  list(file="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/GSE65858_Result.csv",
       prefix="C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/GSE65858",
       escore_color=rgb(148,181,216,255),
       mscore_color=rgb(213,123,112,255),
       emtscore_color=rgb(172,210,199,255))
)
for (dataset in datasets) process_dataset(dataset$file, dataset$prefix, dataset$escore_color, dataset$mscore_color, dataset$emtscore_color)

library(ggpubr)
library(ggplot2)
library(reshape2)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/GSE65858_Result.csv"
data <- read.csv(file_path)

long_data <- melt(data, id.vars = "WDR54", measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
                  variable.name = "SCORE_TYPE", value.name = "SCORE")

custom_palette <- c(
  "E_SCORE" = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE" = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)

gpDensity <- ggscatterhist(
  long_data, x = "WDR54", y = "SCORE",
  color = "SCORE_TYPE", size = 3, alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)

output_file <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/Modified_GSE65858_Result_Plot.png"
ggsave(output_file, plot = gpDensity, width = 9 * 1.2, height = 6)

cat("图形已保存至：", output_file, "\n")

library(ggpubr)
library(ggplot2)
library(reshape2)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/CCLE_Result.csv"
data <- read.csv(file_path)

long_data <- melt(data, id.vars = "WDR54", measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
                  variable.name = "SCORE_TYPE", value.name = "SCORE")

custom_palette <- c(
  "E_SCORE" = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE" = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)

gpDensity <- ggscatterhist(
  long_data, x = "WDR54", y = "SCORE",
  color = "SCORE_TYPE", size = 3, alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)

output_file <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/Modified_CCLE_Result_Plot.png"
ggsave(output_file, plot = gpDensity, width = 9 * 1.2, height = 6)

cat("图形已保存至：", output_file, "\n")

library(ggpubr)
library(ggplot2)
library(reshape2)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/TCGA_Result.csv"
data <- read.csv(file_path)

long_data <- melt(data, id.vars = "WDR54", measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
                  variable.name = "SCORE_TYPE", value.name = "SCORE")

custom_palette <- c(
  "E_SCORE" = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE" = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)

gpDensity <- ggscatterhist(
  long_data, x = "WDR54", y = "SCORE",
  color = "SCORE_TYPE", size = 3, alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)

output_file <- "C:/R/C_Project/F5/Correlation_Heatmap/ssGSEA/results/Modified_TCGA_Result_Plot.png"
ggsave(output_file, plot = gpDensity, width = 9 * 1.2, height = 6)

cat("图形已保存至：", output_file, "\n")

library(dplyr)

gene_list_path <- "C:/R/C_Project/F5/Correlation_Heatmap/gene_list.csv"
ccle_data_path <- "C:/R/C_Project/F5/Correlation_Heatmap/CCLE_HNSCC_ESCC.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/CCLE_Mesenchymal.csv"

gene_list <- read.csv(gene_list_path, header = TRUE, stringsAsFactors = FALSE)

mesenchymal_genes <- gene_list %>%
  filter(Phenotype..marker. == "Mesenchymal") %>%
  select(Symbol) %>%
  pull()

ccle_data <- read.csv(ccle_data_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

selected_genes <- c(mesenchymal_genes, "WDR54")
filtered_data <- ccle_data[, selected_genes, drop = FALSE]

write.csv(filtered_data, file = output_path, row.names = TRUE)

cat("筛选并保存完成，文件保存在:", output_path, "\n")


library(dplyr)

gene_list_path <- "C:/R/C_Project/F5/Correlation_Heatmap/gene_list.csv"
ccle_data_path <- "C:/R/C_Project/F5/Correlation_Heatmap/CCLE_HNSCC_ESCC.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/CCLE_Epithelial.csv"

gene_list <- read.csv(gene_list_path, header = TRUE, stringsAsFactors = FALSE)

epithelial_genes <- gene_list %>%
  filter(Phenotype..marker. == "Epithelial") %>%
  select(Symbol) %>%
  pull()

ccle_data <- read.csv(ccle_data_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

selected_genes <- intersect(c(epithelial_genes, "WDR54"), colnames(ccle_data))
filtered_data <- ccle_data[, selected_genes, drop = FALSE]

write.csv(filtered_data, file = output_path, row.names = TRUE)

cat("筛选并保存完成，文件保存在:", output_path, "\n")


library(dplyr)

gene_list_path <- "C:/R/C_Project/F5/Correlation_Heatmap/gene_list.csv"
tcga_data_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA-HNSC.htseq_counts_annotated.csv"
output_mesenchymal_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA_Mesenchymal.csv"
output_epithelial_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA_Epithelial.csv"

gene_list <- read.csv(gene_list_path, header = TRUE, stringsAsFactors = FALSE)

mesenchymal_genes <- gene_list %>%
  filter(Phenotype..marker. == "Mesenchymal") %>%
  select(Symbol) %>%
  pull()

epithelial_genes <- gene_list %>%
  filter(Phenotype..marker. == "Epithelial") %>%
  select(Symbol) %>%
  pull()

tcga_data <- read.csv(tcga_data_path, header = TRUE, stringsAsFactors = FALSE)

tcga_data_unique <- tcga_data %>%
  distinct(gene_name, .keep_all = TRUE)

selected_mesenchymal_genes <- c(mesenchymal_genes, "WDR54")
mesenchymal_data <- tcga_data_unique %>%
  filter(gene_name %in% selected_mesenchymal_genes)

write.csv(mesenchymal_data, output_mesenchymal_path, row.names = FALSE)

selected_epithelial_genes <- c(epithelial_genes, "WDR54")
epithelial_data <- tcga_data_unique %>%
  filter(gene_name %in% selected_epithelial_genes)

write.csv(epithelial_data, output_epithelial_path, row.names = FALSE)

cat("提取并保存了Mesenchymal和Epithelial基因数据：", output_mesenchymal_path, "和", output_epithelial_path, "\n")


library(dplyr)

gene_list_path <- "C:/R/C_Project/F5/Correlation_Heatmap/gene_list.csv"
gse65858_data_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858.csv"
output_mesenchymal_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858_Mesenchymal.csv"
output_epithelial_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858_Epithelial.csv"

gene_list <- read.csv(gene_list_path, header = TRUE, stringsAsFactors = FALSE)

mesenchymal_genes <- gene_list %>%
  filter(Phenotype..marker. == "Mesenchymal") %>%
  select(Symbol) %>%
  pull()

epithelial_genes <- gene_list %>%
  filter(Phenotype..marker. == "Epithelial") %>%
  select(Symbol) %>%
  pull()

gse65858_data <- read.csv(gse65858_data_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

gse65858_data_unique <- gse65858_data %>%
  distinct(.keep_all = TRUE)

selected_mesenchymal_genes <- c(mesenchymal_genes, "WDR54")
mesenchymal_data <- gse65858_data_unique %>%
  filter(row.names(gse65858_data_unique) %in% selected_mesenchymal_genes)

write.csv(mesenchymal_data, output_mesenchymal_path, row.names = TRUE)

selected_epithelial_genes <- c(epithelial_genes, "WDR54")
epithelial_data <- gse65858_data_unique %>%
  filter(row.names(gse65858_data_unique) %in% selected_epithelial_genes)

write.csv(epithelial_data, output_epithelial_path, row.names = TRUE)

cat("提取并保存了Mesenchymal和Epithelial基因数据：", output_mesenchymal_path, "和", output_epithelial_path, "\n")

normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/CCLE_EM_SCORE.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_CCLE_EM_SCORE.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858_EM_SCORE.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_GSE65858_EM_SCORE.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)
cat("归一化处理完成，并保存到文件：", output_path, "\n")

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA_EM_SCORE.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_TCGA_EM_SCORE.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA_Epithelial.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_TCGA_Epithelial.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
normalize <- function(x) {
  return(2 * (x - min(x)) / (max(x) - min(x)) - 1)
}
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/TCGA_Mesenchymal.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_TCGA_Mesenchymal.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858_Epithelial.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_GSE65858_Epithelial.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

file_path <- "C:/R/C_Project/F5/Correlation_Heatmap/GSE65858_Mesenchymal.csv"
output_path <- "C:/R/C_Project/F5/Correlation_Heatmap/Normalization_GSE65858_Mesenchymal.csv"
data <- read.csv(file_path, header = TRUE, row.names = 1)
data_normalized <- as.data.frame(lapply(data, normalize))
row.names(data_normalized) <- row.names(data)
colnames(data_normalized) <- colnames(data)
write.csv(data_normalized, file = output_path, row.names = TRUE)

library(ComplexHeatmap)
library(circlize)

color_mapping <- colorRamp2(
  c(-2, 0, 2),
  c(rgb(128, 177, 210, maxColorValue = 255),
    rgb(255, 255, 255, maxColorValue = 255),
    rgb(251, 128, 116, maxColorValue = 255))
)

setwd("C:/R/C_Project/F5/Correlation_Heatmap")

file_name <- "Normalization_CCLE_Epithelial.csv"
A <- read.csv(file_name, header = TRUE, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename),
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()

file_name <- "Normalization_CCLE_Mesenchymal.csv"
A <- read.csv(file_name, header = T, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename), 
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()

file_name <- "Normalization_TCGA_Epithelial.csv"
A <- read.csv(file_name, header = T, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename), 
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()

file_name <- "Normalization_TCGA_Mesenchymal.csv"
A <- read.csv(file_name, header = T, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename), 
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()

library(ComplexHeatmap)
library(circlize)

color_mapping <- colorRamp2(
  c(-1, 0, 1), 
  c(
    rgb(128, 177, 210, maxColorValue = 255),  
    rgb(255, 255, 255, maxColorValue = 255),  
    rgb(251, 128, 116, maxColorValue = 255)
  )
)

process_wdr54_heatmap <- function(file_name, output_height, output_width) {
  setwd("C:/R/C_Project/F5/Correlation_Heatmap")
  clinical_data <- read.csv(
    file_name, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1
  )
  wdr54_data <- as.matrix(
    clinical_data[nrow(clinical_data), , drop = FALSE]
  )
  heatmap_plot <- Heatmap(
    wdr54_data, 
    col = color_mapping, 
    show_row_names = FALSE, 
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = NULL, show = FALSE)
  )
  png_filename <- gsub(".csv", "_WDR54_heatmap.png", file_name)
  output_path <- paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename)
  png(
    filename = output_path, 
    width = output_width, 
    height = output_height * 0.5
  )
  draw(heatmap_plot)
  dev.off()
}

process_wdr54_heatmap("Normalization_CCLE_Epithelial.csv", 80, 1000)
process_wdr54_heatmap("Normalization_CCLE_Mesenchymal.csv", 80, 1000)
process_wdr54_heatmap("Normalization_TCGA_Epithelial.csv", 100, 1000)
process_wdr54_heatmap("Normalization_TCGA_Mesenchymal.csv", 100, 1000)

library(ComplexHeatmap)
library(circlize)

color_mapping <- colorRamp2(
  c(-2, 0, 2), 
  c(
    rgb(128, 177, 210, maxColorValue = 255), 
    rgb(255, 255, 255, maxColorValue = 255), 
    rgb(251, 128, 116, maxColorValue = 255)
  )
)

setwd("C:/R/C_Project/F5/Correlation_Heatmap")

file_name <- "Normalization_GSE65858_Epithelial.csv"
A <- read.csv(file_name, header = T, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename), 
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()

file_name <- "Normalization_GSE65858_Mesenchymal.csv"
A <- read.csv(file_name, header = T, row.names = 1)
A <- A[-nrow(A), ]
A <- as.matrix(A)
A[A > 2] <- 2
A[A < -2] <- -2
png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename), 
  width = 900, height = 1200 * 0.8
)
pheatmap(
  A, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
  col = color_mapping, cluster_cols = FALSE
)
dev.off()
library(ComplexHeatmap)
library(circlize)

color_mapping <- colorRamp2(
  c(-1, 0, 1),
  c(
    rgb(128, 177, 210, maxColorValue = 255),
    rgb(255, 255, 255, maxColorValue = 255),
    rgb(251, 128, 116, maxColorValue = 255)
  )
)

process_wdr54_heatmap <- function(file_name, output_height, output_width) {
  setwd("C:/R/C_Project/F5/Correlation_Heatmap")
  clinical_data <- read.csv(
    file_name, header = TRUE, sep = ",",
    stringsAsFactors = FALSE, row.names = 1
  )
  wdr54_data <- as.matrix(clinical_data[nrow(clinical_data), , drop = FALSE])
  heatmap_plot <- Heatmap(
    wdr54_data,
    col = color_mapping,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = NULL, show = FALSE)
  )
  png_filename <- gsub(".csv", "_WDR54_heatmap.png", file_name)
  output_path <- paste0("C:/R/C_Project/F5/Correlation_Heatmap/", png_filename)
  png(
    filename = output_path,
    width = output_width,
    height = output_height * 0.5
  )
  draw(heatmap_plot)
  dev.off()
}

process_wdr54_heatmap("Normalization_GSE65858_Epithelial.csv", 80, 1000)
process_wdr54_heatmap("Normalization_GSE65858_Mesenchymal.csv", 80, 1000)

library(ggpubr)
library(ggplot2)
library(reshape2)

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_CCLE_EM_SCORE.csv")
long_data <- melt(
  data, id.vars = "WDR54",
  measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
  variable.name = "SCORE_TYPE",
  value.name = "SCORE"
)
custom_palette <- c(
  "E_SCORE" = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE" = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)
gpDensity <- ggscatterhist(
  long_data, x = "WDR54", y = "SCORE",
  color = "SCORE_TYPE", size = 3, alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)
ggsave(
  "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Modified_Plot.png",
  plot = gpDensity,
  width = 9 * 1.2, height = 6
)

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_TCGA_EM_SCORE.csv")
long_data <- melt(
  data, id.vars = "WDR54",
  measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
  variable.name = "SCORE_TYPE",
  value.name = "SCORE"
)
custom_palette <- c(
  "E_SCORE" = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE" = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)
gpDensity <- ggscatterhist(
  long_data, x = "WDR54", y = "SCORE",
  color = "SCORE_TYPE", size = 3, alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)
ggsave(
  "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Modified_Plot.png",
  plot = gpDensity,
  width = 9 * 1.2, height = 6
)
library(ggpubr)
library(ggplot2)
library(reshape2)

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_GSE65858_EM_SCORE.csv")
long_data <- melt(data,
                  id.vars = "WDR54",
                  measure.vars = c("E_SCORE", "M_SCORE", "EMT_SCORE"),
                  variable.name = "SCORE_TYPE",
                  value.name = "SCORE")
custom_palette <- c(
  "E_SCORE"   = rgb(148, 181, 216, maxColorValue = 255),
  "M_SCORE"   = rgb(213, 123, 112, maxColorValue = 255),
  "EMT_SCORE" = rgb(172, 210, 199, maxColorValue = 255)
)
gpDensity <- ggscatterhist(
  long_data,
  x = "WDR54",
  y = "SCORE",
  color = "SCORE_TYPE",
  size = 3,
  alpha = 0.6,
  palette = custom_palette,
  margin.params = list(fill = "SCORE_TYPE", color = "black", size = 0.3),
  ggtheme = theme_minimal()
)
ggsave("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Modified_GSE65858_Plot.png",
       plot = gpDensity,
       width = 9 * 1.2,
       height = 6)

calculate_correlation <- function(x_var, y_var, data) {
  pearson_test  <- cor.test(data[[x_var]], data[[y_var]], method = "pearson")
  spearman_test <- cor.test(data[[x_var]], data[[y_var]], method = "spearman")
  data.frame(
    Variable          = y_var,
    Pearson_R         = pearson_test$estimate,
    Pearson_P_value   = pearson_test$p.value,
    Spearman_R        = spearman_test$estimate,
    Spearman_P_value  = spearman_test$p.value
  )
}

plot_and_save <- function(data, x_var, y_var, color, output_file) {
  ggsave(output_file,
         plot = ggscatter(data,
                          x     = x_var,
                          y     = y_var,
                          add   = "reg.line",
                          conf.int = TRUE,
                          color = color,
                          size  = 3,
                          alpha = 0.6,
                          ggtheme = theme_minimal()),
         width  = 9,
         height = 6)
}

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_CCLE_EM_SCORE.csv")
result_e   <- calculate_correlation("WDR54", "E_SCORE", data)
result_m   <- calculate_correlation("WDR54", "M_SCORE", data)
result_emt <- calculate_correlation("WDR54", "EMT_SCORE", data)
final_result <- rbind(result_e, result_m, result_emt)
write.csv(final_result,
          file = "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/CCLE_Correlation_Results.csv",
          row.names = FALSE)
plot_and_save(data, "WDR54", "E_SCORE", rgb(148, 181, 216, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/CCLE_ESCORE.png")
plot_and_save(data, "WDR54", "M_SCORE", rgb(213, 123, 112, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/CCLE_MSCORE.png")
plot_and_save(data, "WDR54", "EMT_SCORE", rgb(172, 210, 199, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/CCLE_EMTSCORE.png")

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_TCGA_EM_SCORE.csv")
result_e   <- calculate_correlation("WDR54", "E_SCORE", data)
result_m   <- calculate_correlation("WDR54", "M_SCORE", data)
result_emt <- calculate_correlation("WDR54", "EMT_SCORE", data)
final_result <- rbind(result_e, result_m, result_emt)
write.csv(final_result,
          file = "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/TCGA_Correlation_Results.csv",
          row.names = FALSE)
plot_and_save(data, "WDR54", "E_SCORE", rgb(148, 181, 216, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/TCGA_ESCORE.png")
plot_and_save(data, "WDR54", "M_SCORE", rgb(213, 123, 112, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/TCGA_MSCORE.png")
plot_and_save(data, "WDR54", "EMT_SCORE", rgb(172, 210, 199, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/TCGA_EMTSCORE.png")

data <- read.csv("C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/Normalization_GSE65858_EM_SCORE.csv")
result_e   <- calculate_correlation("WDR54", "E_SCORE", data)
result_m   <- calculate_correlation("WDR54", "M_SCORE", data)
result_emt <- calculate_correlation("WDR54", "EMT_SCORE", data)
final_result <- rbind(result_e, result_m, result_emt)
write.csv(final_result,
          file = "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/GSE65858_Correlation_Results.csv",
          row.names = FALSE)
plot_and_save(data, "WDR54", "E_SCORE", rgb(148, 181, 216, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/GSE65858_ESCORE.png")
plot_and_save(data, "WDR54", "M_SCORE", rgb(213, 123, 112, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/GSE65858_MSCORE.png")
plot_and_save(data, "WDR54", "EMT_SCORE", rgb(172, 210, 199, maxColorValue = 255),
              "C:/R/C_Project/F5/Correlation_Heatmap/CORRELATION/GSE65858_EMTSCORE.png")



