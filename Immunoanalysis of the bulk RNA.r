### Bulk immune ananlysis
## TIDE
setwd("C:/R/C_Project/F6/TIDE")

data <- read.csv("TCGA_HNSC_mRNA_clean.csv", header = TRUE, sep = ",")

expression_data <- as.matrix(data[, -1])
rownames(expression_data) <- data[, 1]

library(limma)

exprSet <- normalizeBetweenArrays(expression_data)

min_val <- apply(exprSet, 2, min, na.rm = TRUE)
max_val <- apply(exprSet, 2, max, na.rm = TRUE)
normalized_data <- 2 * ((exprSet - min_val) / (max_val - min_val)) - 1

data[, -1] <- normalized_data

output_path_csv_normalized <- "C:/R/C_Project/F6/TIDE/TCGA_HNSC_mRNA_clean_normalized.csv"
write.csv(data, file = output_path_csv_normalized, row.names = FALSE)
print(paste("Min-Max归一化后的CSV文件已保存至:", output_path_csv_normalized))

output_path_tab_normalized <- "C:/R/C_Project/F6/TIDE/TCGA_HNSC_mRNA_clean_normalized.txt"
write.table(data, file = output_path_tab_normalized, sep = "\t", row.names = FALSE, col.names = TRUE)
print(paste("Min-Max归一化后的TAB分隔符文件已保存至:", output_path_tab_normalized))
## 随后在TIDE在线网站上进行运算
setwd("C:/R/C_Project/F6/TIDE")

tcga_tide_data  <- read.csv("TCGA_TIDE.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
tcga_group_data <- read.csv("TCGA_GROUP.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

common_samples     <- intersect(tcga_tide_data[, 1], tcga_group_data[, 1])
tcga_tide_common   <- tcga_tide_data[tcga_tide_data[, 1] %in% common_samples, ]
tcga_group_common  <- tcga_group_data[tcga_group_data[, 1] %in% common_samples, ]
tcga_group_common  <- tcga_group_common[match(tcga_tide_common[, 1], tcga_group_common[, 1]), ]
merged_data        <- cbind(tcga_tide_common, tcga_group_common[, -1])
write.csv(merged_data, file = "WDR54_TIDE.csv", row.names = FALSE)
print("合并后的文件已保存为WDR54_TIDE.csv")

setwd("C:/R/C_Project/F6/TIDE")

data   <- read.csv("WDR54_TIDE.csv", header = TRUE, sep = ",")
wdr54  <- data$WDR54
scores <- data[, 2:12]

results <- data.frame(
  Score       = colnames(scores),
  Pearson_R   = numeric(ncol(scores)),
  Pearson_P   = numeric(ncol(scores)),
  Spearman_R  = numeric(ncol(scores)),
  Spearman_P  = numeric(ncol(scores))
)

for (i in 1:ncol(scores)) {
  score_column          <- scores[, i]
  pearson_test          <- cor.test(wdr54, score_column, method = "pearson")
  results$Pearson_R[i]  <- pearson_test$estimate
  results$Pearson_P[i]  <- pearson_test$p.value
  spearman_test         <- cor.test(wdr54, score_column, method = "spearman")
  results$Spearman_R[i] <- spearman_test$estimate
  results$Spearman_P[i] <- spearman_test$p.value
}

write.csv(results, file = "C:/R/C_Project/F6/TIDE/WDR54_TIDE_correlation_results.csv", row.names = FALSE)
print("相关性计算完成，结果已保存到CSV文件。")

setwd("C:/R/C_Project/F6/TIDE")

data   <- read.csv("WDR54_TIDE.csv", header = TRUE, sep = ",")
group  <- data$Group
scores <- data[, 2:12]

results <- data.frame(
  Score       = colnames(scores),
  T_test_P    = numeric(ncol(scores)),
  Wilcoxon_P  = numeric(ncol(scores))
)

for (i in 1:ncol(scores)) {
  score_column   <- scores[, i]
  group_low      <- score_column[group == "Low"]
  group_high     <- score_column[group == "High"]
  t_test         <- t.test(group_low, group_high)
  results$T_test_P[i]    <- t_test$p.value
  wilcoxon_test         <- wilcox.test(group_low, group_high)
  results$Wilcoxon_P[i]  <- wilcoxon_test$p.value
}

write.csv(results, file = "C:/R/C_Project/F6/TIDE/Group_Significance_Results.csv", row.names = FALSE)
print("显著性差异计算完成，结果已保存到CSV文件。")
library(ggplot2)
library(gghalves)

file_path <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE.csv"
data <- read.csv(file_path)
data$Group <- factor(data$Group, levels = c("High", "Low"))

mycolor_violin <- c("High" = "#D77071", "Low" = "#6888F5")
mycolor_boxplot <- c("High" = "#EF767B", "Low" = "#43A3EF")
mycolor_points <- c("High" = "#C47070", "Low" = "#4370B4")

plot_columns <- colnames(data)[2:12]

plot_violin_box <- function(df, variable, group_var, output_dir) {
    plot_data <- df[, c(variable, group_var)]
    colnames(plot_data) <- c("Value", "Group")
    p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group, color = Group)) +
        geom_point(aes(x = 1.20, y = Value, color = Group),
                   position = position_jitter(width = 0.03),
                   size = 3, shape = 20, alpha = 0.4) +
        geom_boxplot(aes(x = 1.4),
                     position = position_nudge(x = c(0.06, -0.06)),
                     outlier.shape = NA, width = 0.1, alpha = 0.7, color = "black") +
        geom_half_violin(aes(x = 1.6),
                         position = position_nudge(x = -0.05),
                         side = "R", adjust = 1.2, trim = FALSE,
                         color = "black", alpha = 0.8) +
        scale_color_manual(values = mycolor_points) +
        scale_fill_manual(values = mycolor_violin) +
        labs(title = paste("Overlapping", variable, "by Group"),
             x = "Group", y = variable) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5))
    print(p)
    output_file <- file.path(output_dir,
                             paste0(variable, "_Overlap_Violin_Boxplot.png"))
    ggsave(filename = output_file, plot = p, width = 2, height = 6, dpi = 300)
}

output_dir <- "C:/R/C_Project/F6/TIDE/Plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (variable in plot_columns) {
    plot_violin_box(data, variable, "Group", output_dir)
}

library(ggplot2)
library(data.table)

input_file <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE.csv"
output_dir <- "C:/R/C_Project/F6/TIDE/"
data <- fread(input_file)
wdr54_label <- "WDR54"
score_columns <- colnames(data)[2:12]

for (col in score_columns) {
    p <- ggplot(data, aes_string(x = wdr54_label, y = col)) +
        geom_point(size = 2, alpha = 0.6, color = "black") +
        geom_rug(alpha = 0.5, position = "jitter") +
        geom_smooth(method = "lm", se = TRUE,
                    color = "black", linetype = "solid") +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_blank()) +
        labs(x = wdr54_label, y = col)
    output_file <- paste0(output_dir,
                          "Scatter_", col, "_vs_", wdr54_label, ".png")
    ggsave(output_file, plot = p, width = 6, height = 6, dpi = 300)
    cat("图形已保存为:", output_file, "\n")
}

library(data.table)

input_file <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE.csv"
output_file <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE_Selected.csv"
data <- fread(input_file)
selected_columns <- c(1, 2, 9, 10, 13, 14)
data_selected <- data[, ..selected_columns]
normalize_columns <- c(2, 3, 4)
normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}
data_selected[, (normalize_columns) := lapply(.SD, normalize),
              .SDcols = normalize_columns]
fwrite(data_selected, output_file)
cat("归一化后的数据已保存为:", output_file, "\n")

library(ggplot2)
library(data.table)

input_file <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE_Selected.csv"
output_file <- "C:/R/C_Project/F6/TIDE/Integrated_Scatter_Plot.png"
data <- fread(input_file)
wdr54_label <- colnames(data)[5]
group_label <- colnames(data)[6]
score_columns <- colnames(data)[2:4]
data_long <- melt(data,
                  id.vars = c(wdr54_label, group_label),
                  measure.vars = score_columns,
                  variable.name = "Score_Type",
                  value.name = "Score_Value")
mycolor_points <- c("Exclusion" = "#EECA40",
                    "MDSC"      = "#9970AC",
                    "TIDE"      = "#23BAC5")
mycolor_rug <- c("High" = "#D77071", "Low" = "#6888F5")
p <- ggplot(data_long,
            aes_string(x = wdr54_label, y = "Score_Value",
                       color = "Score_Type")) +
    geom_point(size = 2, alpha = 0.6) +
    geom_rug(aes_string(color = group_label), sides = "l",
             alpha = 0.5, position = "jitter") +
    geom_rug(sides = "b", color = "black", alpha = 0.5,
             position = "jitter") +
    geom_smooth(method = "lm", se = TRUE, linetype = "solid") +
    scale_color_manual(values = c(mycolor_points, mycolor_rug)) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_blank()) +
    labs(x = wdr54_label, y = "Score Value")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
cat("整合的散点图已保存为:", output_file, "\n")

library(ggplot2)
library(gghalves)

csv_file_path <- "C:/R/C_Project/F6/TIDE/WDR54_TIDE_Selected_nephelogram.csv"
data <- read.csv(csv_file_path)
colnames(data)

mycolor_violin <- c("TIDE" = "#23BAC5", "Exclusion" = "#EECA40", "MDSC" = "#9970AC")
mycolor_points <- c("TIDE" = "#23BAC5", "Exclusion" = "#EECA40", "MDSC" = "#9970AC")

data$Source <- factor(data$Source, levels = c("TIDE", "Exclusion", "MDSC"))

plot_violin_box <- function(df, variable, group_var, output_dir) {
  plot_data <- df[, c(variable, group_var)]
  colnames(plot_data) <- c("Value", "Source")
  
  p <- ggplot(plot_data, aes(x = Source, y = Value, fill = Source, color = Source)) +
    geom_point(
      aes(x = 1.20, y = Value, color = Source),
      position = position_jitter(width = 0.03),
      size = 0.5,
      shape = 20,
      alpha = 0.4
    ) +
    geom_boxplot(
      aes(x = case_when(
        Source == "TIDE" ~ 1.33,
        Source == "Exclusion" ~ 1.4,
        Source == "MDSC" ~ 1.47
      )),
      outlier.shape = NA,
      width = 0.2,
      alpha = 0.7,
      color = "black"
    ) +
    geom_half_violin(
      aes(x = 1.6),
      position = position_nudge(x = -0.05),
      side = 'R',
      adjust = 1.2,
      trim = FALSE,
      color = "black",
      alpha = 0.8
    ) +
    scale_color_manual(values = mycolor_points) +
    scale_fill_manual(values = mycolor_violin) +
    labs(
      title = "Overlapping Values by Source",
      x = "Source",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  print(p)
  
  output_file <- file.path(output_dir, "WDR54_TIDE_selected.png")
  ggsave(filename = output_file, plot = p, width = 2, height = 6, dpi = 300)
}

output_dir <- "C:/R/C_Project/F6/TIDE"
plot_violin_box(data, "Value", "Source", output_dir)

plot_violin_box <- function(df, variable, group_var, comparisons, output_dir) {
  plot_data <- df[, c(variable, group_var)]
  colnames(plot_data) <- c("Value", "Group")
  
  bg_color <- switch(variable,
    "TIDE" = "#23BAC5",
    "Exclusion" = "#EECA40",
    "MDSC" = "#9970AC",
    "white"
  )
  
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group, color = Group)) +
    geom_violin(alpha = 0.2) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(width = .1, notch = TRUE, alpha = 0.9) +
    geom_quasirandom(aes(color = Group), shape = 16, size = 1, width = 0.1, alpha = 0.4) +
    stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
    labs(x = variable, y = variable) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_color_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
    scale_fill_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      text = element_text(size = 9, color = "black"),
      panel.background = element_rect(fill = alpha(bg_color, 0.1), color = NA)
    )
  
  output_file <- file.path(output_dir, paste0(variable, "_Plot.png"))
  ggsave(filename = output_file, plot = p, width = 4, height = 4, dpi = 300)
}

output_dir <- "C:/R/C_Project/F6/TIDE"
for (variable in columns_to_plot) {
  plot_violin_box(plotdata, variable, "Group", comparisons, output_dir)
}
## TIMER
file_path <- "C:/R/C_Project/F6/TIMER/Infiltration/HNSC_Immune_Infiltration.csv"
data <- read.csv(file_path)

split_info <- strsplit(as.character(data$infiltrates), "_")

data$infiltrates <- sapply(split_info, `[`, 1)
data$database <- sapply(split_info, `[`, 2)

head(data)

write.csv(data, file = "C:/R/C_Project/F6/TIMER/Infiltration/HNSC_Immune_Infiltration_Updated.csv", row.names = FALSE)

cat("文件已成功更新并保存为 HNSC_Immune_Infiltration_Updated.csv。\n")
library(ggplot2)
library(dplyr)

file_path <- "C:/R/C_Project/F6/TIMER/Infiltration/HNSC_Immune_Infiltration_Cleaned.csv"
data <- read.csv(file_path)

data$infiltrates <- factor(data$infiltrates, levels = rev(sort(unique(data$infiltrates))))

data <- data %>%
    mutate(area_scale = case_when(
        p > 0.05 ~ 0.25,
        p <= 0.05 & p > 0.01 ~ 0.25,
        p <= 0.01 & p > 0.001 ~ 0.50,
        p <= 0.001 & p > 0.0001 ~ 0.75,
        p <= 0.0001 ~ 1
    ))

border_color <- "grey"

p <- ggplot(data, aes(x = cancer, y = infiltrates)) +
    geom_tile(color = border_color, fill = NA) +
    geom_tile(aes(fill = rho, width = area_scale, height = area_scale), color = NA) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
    labs(x = "Cancer Type", y = "Immune Infiltrates", fill = "Correlation (Rho)") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, color = "black", size = rel(1.3)),
        axis.text.y.right = element_text(color = "black", size = rel(1.3)),
        axis.title.x = element_text(color = "black", hjust = 0.5),
        axis.title.y.right = element_text(color = "black"),
        axis.ticks.length = unit(0, "pt"),
        axis.text.x.top = element_text(color = "black", angle = 45, hjust = 0),
        axis.title.x.top = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right")

p <- p + geom_text(data = filter(data, p > 0.05), aes(label = "×"), color = border_color, size = 3)

output_path <- "C:/R/C_Project/F6/TIMER/Infiltration/HNSC_Immune_Infiltration_Heatmap.png"
ggsave(output_path, plot = p, width = 10, height = 8)
## TMEsignature
library(dplyr)

input_file <- "C:/R/C_Project/F6/TME_signature/geneset/geneset.csv"
output_dir <- "C:/R/C_Project/F6/TME_signature/geneset"

gene_list <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

write_gmt <- function(gene_list, gene_set_name, output_file) {
  gmt_line <- paste(c(gene_set_name, "NA", gene_list), collapse = "\t")
  writeLines(gmt_line, output_file)
}

for (gene_set_name in colnames(gene_list)) {
  gene_set <- na.omit(gene_list[[gene_set_name]])
  output_file <- file.path(output_dir, paste0(gene_set_name, "_genes.gmt"))
  write_gmt(gene_set, gene_set_name, output_file)
}

cat("所有GMT文件已成功生成并保存到:", output_dir, "\n")

library(clusterProfiler)
library(GSVA)
library(data.table)

expr_file <- "C:/R/C_Project/F6/TME_signature/TCGA.csv"
gmt_dir    <- "C:/R/C_Project/F6/TME_signature/geneset"
output_dir <- "C:/R/C_Project/F6/TME_signature/geneset"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)

expr <- fread(expr_file, header = TRUE, data.table = FALSE, check.names = FALSE)
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

library(stats)

ssgsea_result_file <- "C:/R/C_Project/F6/TME_signature/geneset/ssGSEA_Result.csv"
output_file        <- "C:/R/C_Project/F6/TME_signature/geneset/genes_correlation_results.csv"

ssgsea_result   <- read.csv(ssgsea_result_file, header = TRUE)
data_to_analyze <- ssgsea_result[, 2:7]

calculate_pairwise_correlation <- function(data) {
  num_cols    <- ncol(data)
  cor_results <- data.frame(
    Metric1             = character(),
    Metric2             = character(),
    Spearman_Correlation = numeric(),
    Spearman_p_value     = numeric(),
    Pearson_Correlation  = numeric(),
    Pearson_p_value      = numeric(),
    stringsAsFactors     = FALSE
  )
  for (i in 1:(num_cols - 1)) {
    for (j in (i + 1):num_cols) {
      metric1        <- colnames(data)[i]
      metric2        <- colnames(data)[j]
      spearman_test  <- cor.test(data[, i], data[, j], method = "spearman")
      spearman_corr  <- spearman_test$estimate
      spearman_p     <- spearman_test$p.value
      pearson_test   <- cor.test(data[, i], data[, j], method = "pearson")
      pearson_corr   <- pearson_test$estimate
      pearson_p      <- pearson_test$p.value
      cor_results    <- rbind(cor_results, data.frame(
        Metric1             = metric1,
        Metric2             = metric2,
        Spearman_Correlation = spearman_corr,
        Spearman_p_value     = spearman_p,
        Pearson_Correlation  = pearson_corr,
        Pearson_p_value      = pearson_p,
        stringsAsFactors     = FALSE
      ))
    }
  }
  return(cor_results)
}

correlation_results <- calculate_pairwise_correlation(data_to_analyze)
write.csv(correlation_results, output_file, row.names = FALSE)
cat("相关系数和p值结果已保存至：", output_file, "\n")

library(stats)

tgf_file    <- "C:/R/C_Project/F6/TME_signature/geneset/WDR54_TGF.csv"
output_file <- "C:/R/C_Project/F6/TME_signature/geneset/WDR54_TGF_correlation_results.csv"

tgf_data  <- read.csv(tgf_file, header = TRUE)
WDR54     <- tgf_data[, "WDR54"]
tgf_genes <- tgf_data[, 2:12]

calculate_correlation <- function() {
  cor_results <- data.frame(
    Gene                 = colnames(tgf_genes),
    Spearman_Correlation = numeric(length(colnames(tgf_genes))),
    Spearman_p_value      = numeric(length(colnames(tgf_genes))),
    Pearson_Correlation  = numeric(length(colnames(tgf_genes))),
    Pearson_p_value      = numeric(length(colnames(tgf_genes)))
  )
  for (i in 1:ncol(tgf_genes)) {
    spearman_test <- cor.test(WDR54, tgf_genes[, i], method = "spearman")
    cor_results$Spearman_Correlation[i] <- spearman_test$estimate
    cor_results$Spearman_p_value[i]     <- spearman_test$p.value
    pearson_test                         <- cor.test(WDR54, tgf_genes[, i], method = "pearson")
    cor_results$Pearson_Correlation[i]  <- pearson_test$estimate
    cor_results$Pearson_p_value[i]      <- pearson_test$p.value
  }
  return(cor_results)
}

correlation_results <- calculate_correlation()
write.csv(correlation_results, output_file, row.names = FALSE)
cat("相关系数和p值结果已保存至：", output_file, "\n")

library(circlize)

data             <- read.csv("C:/R/C_Project/F6/TME_signature/geneset/circolus.csv")
significant_data <- data[data$Spearman_p_value <= 0.05, ]

sector_colors <- c(
  "INF_gamma_Response"             = "#F07673",
  "Wound_Healing"                  = "#6888F5",
  "Overall_Lymphocyte_Infiltration"= "#9970AC",
  "TGF_beta_Response"              = "#23BAC5",
  "Macrophages_Monocytes_Regulation" = "#7AD151",
  "WDR54"                          = "#F1D756"
)

opacity_values <- ifelse(significant_data$Spearman_p_value == 0, 0.90,
                   ifelse(significant_data$Spearman_p_value < 1e-4, 0.80,
                   ifelse(significant_data$Spearman_p_value < 1e-3, 0.60,
                   ifelse(significant_data$Spearman_p_value < 1e-2, 0.40,
                   ifelse(significant_data$Spearman_p_value < 0.05, 0.20, 0)))))

sector_order <- c(
  "INF_gamma_Response",
  "Wound_Healing",
  "Overall_Lymphocyte_Infiltration",
  "Macrophages_Monocytes_Regulation",
  "WDR54",
  "TGF_beta_Response"
)

png("C:/R/C_Project/F6/TME_signature/geneset/correlation_chord_diagram_with_track.png",
    width = 800, height = 800)

chordDiagram(
  significant_data[, c("Metric1", "Metric2", "Spearman_Correlation")],
  transparency    = 1 - opacity_values,
  grid.col        = sector_colors,
  col             = "#FFE0C1",
  annotationTrack = "grid",
  order           = sector_order,
  preAllocateTracks = list(track.height = 0.2),
  link.lwd        = abs(significant_data$Spearman_Correlation)
)

circos.trackPlotRegion(
  track.index  = 1,
  track.height = 1,
  panel.fun    = function(x, y) {
    circos.axis(labels.cex = 1.2, lwd = 2, major.tick.length = 0.2)
  },
  bg.border    = NA
)

dev.off()
circos.clear()

library(ggplot2)

setwd("C:/R/C_Project/F6/TME_signature/geneset")

bar_data <- read.csv("Bar.csv", header = TRUE)
bar_data$Gene        <- factor(bar_data$Gene, levels = bar_data$Gene)
bar_data$log_p_value <- log10(bar_data$Pearson_p_value)

ggplot(bar_data, aes(x = Gene, y = Pearson_Correlation, fill = log_p_value)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_gradient(
    low    = "#6888F5",
    high   = "#F07673",
    limits = c(log10(1e-40), log10(1e-2)),
    breaks = log10(c(1e-40, 1e-20, 1e-10, 1e-5, 1e-2)),
    labels = c("1e-40", "1e-20", "1e-10", "1e-5", "1e-2")
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x       = element_text(angle = 90, hjust = 1),
    panel.background  = element_rect(fill = "white", color = "white"),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  ) +
  labs(
    title = "Pearson Correlation by Gene",
    x     = "Gene",
    y     = "Pearson Correlation",
    fill  = "Log10 P-value"
  )

ggsave("Bar_plot.png", width = 8, height = 6, dpi = 300)
