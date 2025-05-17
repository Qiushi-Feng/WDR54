### ssGSEA
library(dplyr)

input_file <- "C:/R/C_Project/F6/ssGSEA/gene_list.csv"
output_dir  <- "C:/R/C_Project/F6/ssGSEA/geneset"

gene_list <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
colnames(gene_list) <- c("Symbol", "GeneSet")

write_gmt <- function(gene_list, gene_set_name, output_file) {
  gmt_line <- paste(c(gene_set_name, "NA", gene_list), collapse = "\t")
  writeLines(gmt_line, output_file)
}

gene_sets <- split(gene_list$Symbol, gene_list$GeneSet)

for (gene_set_name in names(gene_sets)) {
  output_file <- file.path(output_dir, paste0(gene_set_name, "_genes.gmt"))
  write_gmt(gene_sets[[gene_set_name]], gene_set_name, output_file)
}

cat("所有GMT文件已成功生成并保存到:", output_dir, "\n")

library(clusterProfiler)
library(GSVA)
library(data.table)

expr_file  <- "C:/R/C_Project/F6/ssGSEA/TCGA.csv"
gmt_dir    <- "C:/R/C_Project/F6/ssGSEA/geneset"
output_dir <- "C:/R/C_Project/F6/ssGSEA/results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

gmt_files   <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
expr        <- fread(expr_file, header = TRUE, data.table = FALSE, check.names = FALSE)
rownames(expr) <- expr[, 1]
expr[, -1]     <- apply(expr[, -1], 2, function(x) as.numeric(as.character(x)))
expr_matrix    <- as.matrix(expr[, -1])

run_ssgsea_for_gmt <- function(gmt_file) {
  genesets       <- clusterProfiler::read.gmt(gmt_file)
  genesets4gsva  <- split(genesets$gene, genesets$term)
  ssgsea_param   <- ssgseaParam(expr_matrix, genesets4gsva, alpha = 0.25, normalize = TRUE)
  ssgsea_res     <- gsva(ssgsea_param)
  gene_set_name  <- gsub(".gmt", "", basename(gmt_file))
  ssgsea_df      <- as.data.frame(ssgsea_res)
  output_file    <- file.path(output_dir, paste0(gene_set_name, "_ESCORE.csv"))
  write.csv(ssgsea_df, output_file, row.names = TRUE)
  cat("结果已保存至：", output_file, "\n")
}

for (gmt_file in gmt_files) {
  run_ssgsea_for_gmt(gmt_file)
}

cat("所有结果已生成并保存至：", output_dir, "\n")
input_file  <- "C:/R/C_Project/F6/ssGSEA/results/ssGSEA_Result.csv"
output_file <- "C:/R/C_Project/F6/ssGSEA/results/ssGSEA_Result_cleaned.csv"

data <- read.csv(input_file, row.names = 1)

data_to_clean <- data[2:34, ]

detect_outliers <- function(row) {
  Q1 <- quantile(row, 0.25, na.rm = TRUE)
  Q3 <- quantile(row, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  (row < lower_bound) | (row > upper_bound)
}

outlier_columns <- c()

for (i in 1:nrow(data_to_clean)) {
  outliers <- detect_outliers(data_to_clean[i, ])
  column_names_with_outliers <- colnames(data_to_clean)[outliers]
  outlier_columns <- unique(c(outlier_columns, column_names_with_outliers))
}

cleaned_data <- data[, !(colnames(data) %in% outlier_columns)]

write.csv(cleaned_data, output_file, row.names = TRUE)

input_file  <- "C:/R/C_Project/F6/ssGSEA/results/ssGSEA_Result_cleaned.csv"
output_file <- "C:/R/C_Project/F6/ssGSEA/results/correlation_results.csv"

data <- read.csv(input_file, row.names = 1)

WDR54_data      <- as.numeric(data[nrow(data), ])
data_to_analyze <- data[1:(nrow(data) - 1), ]

results <- data.frame(
  Row                 = rownames(data_to_analyze),
  Pearson_Correlation = NA,
  Pearson_P_Value     = NA,
  Spearman_Correlation= NA,
  Spearman_P_Value    = NA
)

for (i in 1:nrow(data_to_analyze)) {
  current_row <- as.numeric(data_to_analyze[i, ])
  pearson_test <- cor.test(WDR54_data, current_row, method = "pearson", use = "complete.obs")
  results$Pearson_Correlation[i] <- pearson_test$estimate
  results$Pearson_P_Value[i]       <- pearson_test$p.value
  spearman_test <- cor.test(WDR54_data, current_row, method = "spearman", use = "complete.obs")
  results$Spearman_Correlation[i] <- spearman_test$estimate
  results$Spearman_P_Value[i]     <- spearman_test$p.value
}

write.csv(results, output_file, row.names = FALSE)
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

setwd("C:/R/C_Project/F6/ssGSEA/results")

file_name <- "ssGSEA_Result_cleaned.csv"
A <- read.csv(file_name, header = TRUE, row.names = 1)
A <- as.matrix(A)
A <- A[-nrow(A), ]
A[A > 2]  <- 2
A[A < -2] <- -2

png_filename <- gsub(".csv", "_heatmap.png", file_name)
png(
  filename = paste0("C:/R/C_Project/F6/ssGSEA/results/", png_filename),
  width    = 900,
  height   = 1200 * 0.8
)

pheatmap(
  A,
  scale          = "row",
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  col            = color_mapping,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE
)

dev.off()

cat(
  "热图已生成并保存到:",
  paste0("C:/R/C_Project/F6/ssGSEA/results/", png_filename),
  "\n"
)

library(tidyverse)

file_path <- "C:/R/C_Project/F6/ssGSEA/results/ssGSEA_Result_cleaned.csv"
data      <- read.csv(file_path, header = TRUE, row.names = 1)

wdr54_expr        <- as.numeric(data["WDR54", ])
wdr54_expr        <- wdr54_expr[is.finite(wdr54_expr)]
wdr54_expr_sorted <- sort(wdr54_expr)
mid_index         <- length(wdr54_expr_sorted) %/% 2

colors <- c(
  rep("#6888F5", mid_index),
  rep("#D77071", length(wdr54_expr_sorted) - mid_index)
)

output_file <- "C:/R/C_Project/F6/ssGSEA/results/WDR54_expression_curve.png"
png(filename = output_file, width = 6000, height = 1000)

plot(
  wdr54_expr_sorted,
  type  = "n",
  xlab  = "",
  ylab  = "",
  main  = "",
  axes  = FALSE
)
lines(
  1:mid_index,
  wdr54_expr_sorted[1:mid_index],
  col = "#6888F5",
  lwd = 10
)
lines(
  (mid_index + 1):length(wdr54_expr_sorted),
  wdr54_expr_sorted[(mid_index + 1):length(wdr54_expr_sorted)],
  col = "#D77071",
  lwd = 10
)

axis(1, lwd = 4, cex.axis = 4.0)
axis(2, lwd = 4, cex.axis = 4.0)

dev.off()

cat("WDR54表达曲线图已生成并保存到:", output_file, "\n")

input_dir   <- "C:/R/C_Project/F6/ssGSEA/geneset/"
output_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/geneset/Angiogenesis_Fibroblast.gmt"
geneset_files <- c(
  "Angiogenesis",
  "Cancer-associated fibroblasts",
  "Endothelium",
  "Granulocyte traffic",
  "Matrix",
  "Matrix remodeling",
  "Neutrophil signature",
  "Protumor cytokines"
)
all_genesets <- c()
for (file_name in geneset_files) {
  file_path        <- paste0(input_dir, file_name, "_genes.gmt")
  geneset_content  <- readLines(file_path, warn = FALSE)
  all_genesets     <- c(all_genesets, geneset_content)
}
writeLines(all_genesets, output_file)
cat("基因集文件已整合并保存为:", output_file, "\n")

input_dir   <- "C:/R/C_Project/F6/ssGSEA/geneset/"
output_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/geneset/Pro-tumor_Immune_Infiltrate.gmt"
geneset_files <- c(
  "Immune Suppression by Myeloid Cells",
  "M1 signature",
  "Macrophage and DC traffic",
  "Myeloid cells traffic",
  "Th2 signature",
  "Treg",
  "Treg and Th2 traffic",
  "Tumor-associated Macrophages"
)
all_genesets <- c()
for (file_name in geneset_files) {
  file_path        <- paste0(input_dir, file_name, "_genes.gmt")
  geneset_content  <- readLines(file_path, warn = FALSE)
  all_genesets     <- c(all_genesets, geneset_content)
}
writeLines(all_genesets, output_file)
cat("基因集文件已整合并保存为:", output_file, "\n")

input_dir   <- "C:/R/C_Project/F6/ssGSEA/geneset/"
output_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/geneset/EMT_Proliferation.gmt"
geneset_files <- c("EMT signature", "Tumor proliferation rate")
all_genesets <- c()
for (file_name in geneset_files) {
  file_path        <- paste0(input_dir, file_name, "_genes.gmt")
  geneset_content  <- readLines(file_path, warn = FALSE)
  all_genesets     <- c(all_genesets, geneset_content)
}
writeLines(all_genesets, output_file)
cat("基因集文件已整合并保存为:", output_file, "\n")
library(clusterProfiler)
library(GSVA)
library(data.table)

expr_file <- "C:/R/C_Project/F6/ssGSEA/TCGA.csv"
gmt_dir   <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/geneset"
output_dir<- "C:/R/C_Project/F6/ssGSEA/scatter_diagram"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)

expr <- fread(expr_file, header = TRUE, data.table = FALSE, check.names = FALSE)
rownames(expr) <- expr[,1]
expr[,-1] <- apply(expr[,-1], 2, function(x) as.numeric(as.character(x)))
expr_matrix <- as.matrix(expr[,-1])

run_ssgsea_for_gmt <- function(gmt_file) {
  genesets       <- clusterProfiler::read.gmt(gmt_file)
  genesets4gsva  <- split(genesets$gene, genesets$term)
  ssgsea_param   <- ssgseaParam(expr_matrix, genesets4gsva, alpha = 0.25, normalize = TRUE)
  ssgsea_res     <- gsva(ssgsea_param)
  gene_set_name  <- gsub(".gmt", "", basename(gmt_file))
  ssgsea_df      <- as.data.frame(ssgsea_res)
  output_file    <- file.path(output_dir, paste0(gene_set_name, "_ESCORE.csv"))
  write.csv(ssgsea_df, output_file, row.names = TRUE)
  cat("结果已保存至：", output_file, "\n")
}

for (gmt_file in gmt_files) run_ssgsea_for_gmt(gmt_file)

cat("所有结果已生成并保存至：", output_dir, "\n")

library(data.table)

input_file  <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/ssGSEA_Result.csv"
output_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/correlation_results.csv"

data <- fread(input_file)
wdr54_expr <- data[[6]]

results <- data.frame(
  Method         = character(),
  Score_Column   = character(),
  Correlation    = numeric(),
  P_value        = numeric(),
  stringsAsFactors = FALSE
)

score_columns <- names(data)[2:5]

for (col in score_columns) {
  score_data     <- data[[col]]
  pearson_test   <- cor.test(score_data, wdr54_expr, method = "pearson")
  results        <- rbind(results, data.frame(
    Method       = "Pearson",
    Score_Column = col,
    Correlation  = pearson_test$estimate,
    P_value      = pearson_test$p.value
  ))
  spearman_test <- cor.test(score_data, wdr54_expr, method = "spearman")
  results       <- rbind(results, data.frame(
    Method       = "Spearman",
    Score_Column = col,
    Correlation  = spearman_test$estimate,
    P_value      = spearman_test$p.value
  ))
}

write.csv(results, output_file, row.names = FALSE)
cat("相关性分析已完成，结果保存到:", output_file, "\n")

library(data.table)

input_file  <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/ssGSEA_Result.csv"
output_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/ssGSEA_Result_normalized.csv"

data <- fread(input_file)

normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

data[,2:6] <- lapply(data[,2:6, with = FALSE], normalize)

fwrite(data, output_file)
cat("数据已归一化并保存至:", output_file, "\n")

library(ggplot2)
library(data.table)

input_file <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/ssGSEA_Result_normalized.csv"
output_dir <- "C:/R/C_Project/F6/ssGSEA/scatter_diagram/"

data <- fread(input_file)
wdr54_label <- "WDR54"

score_columns <- c(
  "Pro_tumor_Immune_Infiltrate",
  "Angiogenesis_Fibroblast",
  "Anti_tumor_Immune_Infiltrate",
  "EMT_Proliferation"
)
column_colors <- c(
  "Pro_tumor_Immune_Infiltrate" = "#F1D756",
  "Angiogenesis_Fibroblast"     = "#6888F5",
  "Anti_tumor_Immune_Infiltrate"= "#9970AC",
  "EMT_Proliferation"           = "#F07673"
)

for (col in score_columns) {
  p <- ggplot(data, aes_string(x = wdr54_label, y = col)) +
    geom_point(size = 2, alpha = 0.6, color = column_colors[col]) +
    geom_rug(alpha = 0.5, position = "jitter") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_blank()) +
    labs(x = wdr54_label, y = col)
  output_file <- paste0(output_dir, "Scatter_", col, "_vs_", wdr54_label, ".png")
  ggsave(output_file, plot = p, width = 8 * 2/3, height = 6 * 2/3, dpi = 300)
  cat("图形已保存为:", output_file, "\n")
}
library(igraph)
library(fpc)
library(cluster)

setwd("C:/R/C_Project/F6/Louvain_cluster")

ssGSEA_data <- read.csv("ssGSEA_Result.csv", row.names = 1, check.names = FALSE)

ssGSEA_data_std <- scale(ssGSEA_data, center = TRUE, scale = TRUE)
write.csv(ssGSEA_data_std, "standardized_ssGSEA_Result.csv", row.names = TRUE)

pearson_corr <- cor(t(ssGSEA_data_std), method = "pearson")
write.csv(pearson_corr, "pearson_correlation_matrix.csv", row.names = TRUE)

threshold <- 0.45
adj_matrix <- pearson_corr
adj_matrix[adj_matrix < threshold] <- 0
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
write.csv(adj_matrix, "adjacency_matrix.csv", row.names = TRUE)

louvain_communities <- cluster_louvain(g)
membership_vector <- membership(louvain_communities)
write.csv(membership_vector, "louvain_membership.csv", row.names = TRUE)

dist_matrix <- as.dist(1 - pearson_corr)
sil_scores <- silhouette(membership_vector, dist_matrix)
write.csv(sil_scores[, 1:3], "silhouette_scores.csv", row.names = TRUE)
avg_sil_width <- mean(sil_scores[, 3])

layout_fr <- layout_with_fr(g, niter = 1000, area = vcount(g)^2.5, repulserad = vcount(g)^3)

png("louvain_clusters.png", width = 2400, height = 1800)
plot(louvain_communities, g, layout = layout_fr, vertex.size = 5, vertex.label = NA, main = "Louvain Clustering")
dev.off()

png("silhouette_plot.png", width = 2400, height = 1800)
plot(sil_scores, main = paste("Silhouette Plot (Avg width:", round(avg_sil_width, 2), ")"))
dev.off()

write.graph(g, file = "graph_for_cytoscape.graphml", format = "graphml")

library(pheatmap)

louvain_data <- read.csv("C:/R/C_Project/F6/Louvain_cluster/louvain_membership.csv")
ssgsea_data  <- read.csv("C:/R/C_Project/F6/Louvain_cluster/ssGSEA_Result.csv")

samples <- louvain_data[, 1]
groups  <- louvain_data[, 2]
grouped_samples <- split(samples, groups)

combined_data <- data.frame()
col_counter <- 1

normalize_row <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  if (max_val != min_val) {
    return(2 * (x - min_val) / (max_val - min_val) - 1)
  } else {
    return(rep(0, length(x)))
  }
}
for (group in unique(groups)) {
  group_samples <- grouped_samples[[as.character(group)]]
  group_data <- ssgsea_data[ssgsea_data[, 1] %in% group_samples, ]
  group_data <- group_data[, -ncol(group_data)]
  transposed_data <- t(group_data[, -1])
  transposed_data <- as.data.frame(transposed_data)
  colnames(transposed_data) <- group_data[, 1]
  normalized_data <- t(apply(transposed_data, 1, normalize_row))
  normalized_data <- as.data.frame(normalized_data)
  colnames(normalized_data) <- group_data[, 1]
  if (nrow(combined_data) == 0) {
    combined_data <- normalized_data
  } else {
    for (i in 1:20) {
      combined_data[[sprintf("%03d", col_counter)]] <- 0
      col_counter <- col_counter + 1
    }
    combined_data <- cbind(combined_data, normalized_data)
  }
}

write.csv(combined_data, "C:/R/C_Project/F6/Louvain_cluster/combined_heatmap.csv", row.names = TRUE)

heatmap_matrix <- as.matrix(combined_data)

png("C:/R/C_Project/F6/Louvain_cluster/combined_heatmap.png", width = 600, height = 800)
pheatmap(heatmap_matrix,
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         show_rownames  = TRUE,
         show_colnames  = TRUE)
dev.off()

louvain_data     <- read.csv("C:/R/C_Project/F6/Louvain_cluster/louvain_membership.csv")
expression_data  <- read.csv("C:/R/C_Project/F6/Louvain_cluster/TCGA_HNSC_mRNA_clean.csv", row.names = 1)

samples          <- louvain_data[, 1]
groups           <- louvain_data[, 2]
grouped_samples  <- split(samples, groups)

WDR54_expression <- expression_data["WDR54", , drop = FALSE]

combined_expression <- data.frame(
  Sample            = character(),
  WDR54_Expression  = numeric(),
  Group             = character(),
  stringsAsFactors  = FALSE
)

for (group in names(grouped_samples)) {
  group_samples    <- grouped_samples[[group]]
  valid_samples    <- group_samples[group_samples %in% colnames(WDR54_expression)]
  group_expression <- as.numeric(WDR54_expression[, valid_samples])
  temp_df          <- data.frame(
    Sample           = valid_samples,
    WDR54_Expression = group_expression,
    Group            = group,
    stringsAsFactors = FALSE
  )
  combined_expression <- rbind(combined_expression, temp_df)
}

output_file <- "C:/R/C_Project/F6/Louvain_cluster/WDR54_expression_combined.csv"
write.csv(combined_expression, output_file, row.names = FALSE)
cat("WDR54表达水平已提取并保存为一个CSV文件，包含样本名称、WDR54表达和样本分组。\n")

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

rm(list = ls())
options(stringsAsFactors = FALSE)

file_path <- "C:/R/C_Project/F6/Louvain_cluster/WDR54_expression_combined.csv"
plotdata  <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

head(plotdata)

plotdata     <- plotdata %>% mutate(Group = factor(Group, levels = unique(Group)))
comparisons  <- combn(levels(plotdata$Group), 2, simplify = FALSE)

p_values     <- compare_means(
  WDR54_Expression ~ Group,
  data        = plotdata,
  comparisons = comparisons,
  method      = "t.test"
)

p_values_path <- "C:/R/C_Project/F6/Louvain_cluster/WDR54_expression_P_values.csv"
write.csv(p_values, file = p_values_path, row.names = FALSE)

pl <- ggplot(plotdata, aes(x = Group, y = WDR54_Expression)) +
  geom_violin(aes(fill = Group), alpha = 0.1, color = NA) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = .1, notch = TRUE) +
  geom_quasirandom(
    aes(group = Group),
    shape = 16, size = 0.5, width = 0.1,
    color = c("#6888F5", "#D77071", "#549F9A", "#FD763F")[as.numeric(plotdata$Group)]
  ) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  labs(x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#6888F5", "#D77071", "#549F9A", "#FD763F")) +
  scale_fill_manual(values = c("#6888F5", "#D77071", "#549F9A", "#FD763F")) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F6/Louvain_cluster/WDR54_expression_Plot.png"
ggsave(filename = output_path, plot = pl, width = 5.33, height = 4, units = "in", dpi = 300)
