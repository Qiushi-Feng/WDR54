## Bulk preprocession and debatch
file_path <- "C:/R/C_Project/GEO_output/GSE85446/GSE85446_clean.csv"
exp1 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp1)

file_path <- "C:/R/C_Project/GEO_output/GSE41613/GSE41613_clean.csv"
exp2 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp2)

file_path <- "C:/R/C_Project/GEO_output/GSE65858/GSE65858_clean.csv"
exp3 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp3)

file_path <- "C:/R/C_Project/GEO_output/GSE75538/GSE75538_clean.csv"
exp4 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp4)

file_path <- "C:/R/C_Project/GEO_output/GSE30784/GSE30784_clean.csv"
exp5 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp5)

file_path <- "C:/R/C_Project/TCGA_output/mRNA/TCGA_HNSC_mRNA_clean.csv"
exp6 <- read.csv(file_path, header = TRUE, row.names = 1)
dim(exp6)

colnames(exp1) <- paste0("GSE85446", 1:66)
colnames(exp2) <- paste0("GSE41613", 1:97)
colnames(exp3) <- paste0("GSE65858", 1:270)
colnames(exp4) <- paste0("GSE75538", 1:28)
colnames(exp5) <- paste0("GSE30784", 1:229)
colnames(exp6) <- paste0("TCGA_HNSC", 1:545)

all_genes <- Reduce(union, list(
  rownames(exp1), rownames(exp2), rownames(exp3),
  rownames(exp4), rownames(exp5), rownames(exp6)
))

exp1_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp1),
  dimnames = list(all_genes, colnames(exp1))
)
exp2_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp2),
  dimnames = list(all_genes, colnames(exp2))
)
exp3_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp3),
  dimnames = list(all_genes, colnames(exp3))
)
exp4_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp4),
  dimnames = list(all_genes, colnames(exp4))
)
exp5_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp5),
  dimnames = list(all_genes, colnames(exp5))
)
exp6_full <- matrix(0,
  nrow = length(all_genes), ncol = ncol(exp6),
  dimnames = list(all_genes, colnames(exp6))
)

exp1_full[rownames(exp1), ] <- as.matrix(exp1)
exp2_full[rownames(exp2), ] <- as.matrix(exp2)
exp3_full[rownames(exp3), ] <- as.matrix(exp3)
exp4_full[rownames(exp4), ] <- as.matrix(exp4)
exp5_full[rownames(exp5), ] <- as.matrix(exp5)
exp6_full[rownames(exp6), ] <- as.matrix(exp6)

exp_all <- cbind(
  exp1_full, exp2_full, exp3_full,
  exp4_full, exp5_full, exp6_full
)

zero_counts <- rowSums(exp_all == 0)
threshold <- 0.75 * ncol(exp_all)
genes_to_keep <- zero_counts < threshold
exp_all <- exp_all[genes_to_keep, ]

output_path_filtered <- "C:/R/C_Project/debatch_data/clean_mRNA_50%union.csv"
write.csv(exp_all, file = output_path_filtered, row.names = TRUE)

dim(exp_all)
head(exp_all)

group_list <- data.frame(
  sample = colnames(exp_all), 
  dataset = c(rep("GSE30784", ncol(exp1)), rep("GSE41613", ncol(exp2)), rep("GSE65858", ncol(exp3)), 
              rep("GSE75538", ncol(exp4)), rep("GSE85446", ncol(exp5)), rep("TCGA_HNSC", ncol(exp6)))
)
rownames(group_list) <- group_list$sample

library(umap)

umap_result <- umap(t(exp_all))

umap_data <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  dataset = group_list$dataset
)

library(ggplot2)
umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Gene Expression Data", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

output_path <- "C:/R/C_Project/debatch_data/UMAP_plot_before_debatch.png"
ggsave(output_path, plot = umap_plot, width = 8, height = 6)

cat("UMAP plot saved to:", output_path, "\n")

pca_result <- prcomp(t(exp_all), scale. = TRUE)

pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  dataset = group_list$dataset
)

library(ggplot2)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "PCA of Gene Expression Data", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

output_path <- "C:/R/C_Project/debatch_data/PCA_plot_before_debatch.png"
ggsave(output_path, plot = pca_plot, width = 8, height = 6)

cat("PCA plot saved to:", output_path, "\n")

library(tidyverse)

exp_all_long <- as.data.frame(exp_all) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "expression", -gene)

exp_all_long$batch <- case_when(
  grepl("^GSE85446", exp_all_long$sample) ~ "GSE85446",
  grepl("^GSE41613", exp_all_long$sample) ~ "GSE41613",
  grepl("^GSE65858", exp_all_long$sample) ~ "GSE65858",
  grepl("^GSE75538", exp_all_long$sample) ~ "GSE75538",
  grepl("^GSE30784", exp_all_long$sample) ~ "GSE30784",
  grepl("^TCGA_HNSC", exp_all_long$sample) ~ "TCGA_HNSC"
)
dataset_colors <- c(
  "GSE30784" = rgb(248, 118, 109, maxColorValue = 255),
  "GSE41613" = rgb(183, 159, 0, maxColorValue = 255),
  "GSE65858" = rgb(0, 186, 56, maxColorValue = 255),
  "GSE75538" = rgb(0, 191, 196, maxColorValue = 255),
  "GSE85446" = rgb(97, 156, 255, maxColorValue = 255),
  "TCGA_HNSC" = rgb(245, 100, 227, maxColorValue = 255)
)
exp_all_long$batch <- factor(exp_all_long$batch, levels = names(dataset_colors))

p <- ggplot(exp_all_long, aes(x = sample, y = expression, fill = batch, colour = batch)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  scale_colour_manual(values = dataset_colors) +
  theme_bw() +
  labs(title = "Boxplot of Gene Expression by Sample",
       x = "Sample",
       y = "Gene Expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

file_path <- "C:/R/C_Project/debatch_data/boxplot_before_debatch.png"

ggsave(filename = file_path, plot = p, width = 15, height = 10)

cat("Boxplot saved to:", file_path, "\n")

library(tidyverse)
library(ggdendro)

data_to_cluster <- exp_all

tmp_data <- dist(as.data.frame(t(data_to_cluster)))

clustering <- hclust(tmp_data, method = "average")

dendro_data <- dendro_data(clustering)

dendro_data$labels$batch <- factor(group_list$dataset[match(dendro_data$labels$label, group_list$sample)])

dataset_colors <- c(
  "GSE30784"   = rgb(248, 118, 109, maxColorValue = 255),
  "GSE41613"   = rgb(183, 159,   0, maxColorValue = 255),
  "GSE65858"   = rgb(  0, 186,  56, maxColorValue = 255),
  "GSE75538"   = rgb(  0, 191, 196, maxColorValue = 255),
  "GSE85446"   = rgb( 97, 156, 255, maxColorValue = 255),
  "TCGA_HNSC"  = rgb(245, 100, 227, maxColorValue = 255)
)

p <- ggplot() +
  geom_segment(data = dendro_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dendro_data$labels,
             aes(x = x, y = y, color = batch), size = 3) +
  scale_color_manual(values = dataset_colors) +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.margin      = margin(1, 1, 1, 1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Hierarchical Clustering of Samples") +
  coord_cartesian(ylim = c(0, max(dendro_data$segments$y) * 0.6))



# 设置文件保存路径
file_path_clustering <- "C:/R/C_Project/debatch_data/hierarchical_clustering_before_debatch.png"

# 保存聚类树图形到本地文件
ggsave(filename = file_path_clustering, plot = p, width = 15, height = 10)

# 打印保存图像的路径
cat("Hierarchical clustering plot saved to:", file_path_clustering, "\n")

# 方案一使用SVA包去除批次效应
library(sva)
exp_all_combat <- ComBat(exp_all, batch = group_list$dataset) # batch为批次信息


# 保存去除批次效应后的数据到CSV文件
output_path_debatch <- "C:/R/C_Project/debatch_data/debatch_mRNA_50%union.csv"
write.csv(exp_all_combat, file = output_path_debatch, row.names = TRUE)

# 打印保存去批次效应数据的路径
cat("Debatch data saved to:", output_path_debatch, "\n")


# 设置分组信息
group_list <- data.frame(
  sample = colnames(exp_all_combat), 
  dataset = c(rep("GSE30784", ncol(exp1)), rep("GSE41613", ncol(exp2)), rep("GSE65858", ncol(exp3)), 
              rep("GSE75538", ncol(exp4)), rep("GSE85446", ncol(exp5)), rep("TCGA_HNSC", ncol(exp6)))
)
rownames(group_list) <- group_list$sample


## UMAP分析（推荐）
library(umap)

# 进行 UMAP 分析
umap_result_combat <- umap(t(exp_all_combat))

# 提取 UMAP 结果
umap_data_combat <- data.frame(
  UMAP1 = umap_result_combat$layout[, 1],
  UMAP2 = umap_result_combat$layout[, 2],
  dataset = group_list$dataset
)

# 可视化 UMAP 结果
umap_plot_combat <- ggplot(umap_data_combat, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Gene Expression Data After Batch Effect Removal", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# 保存 UMAP 图像
output_path_combat <- "C:/R/C_Project/debatch_data/UMAP_plot_after_debatch.png"
ggsave(output_path_combat, plot = umap_plot_combat, width = 8, height = 6)

# 打印保存图像的路径
cat("UMAP plot after batch effect removal saved to:", output_path_combat, "\n")


## PCA主成分分析

# 进行 PCA 分析
pca_result <- prcomp(t(exp_all_combat), scale. = TRUE)

# 提取前两个主成分
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  dataset = group_list$dataset
)
file_path_clustering <- "C:/R/C_Project/debatch_data/hierarchical_clustering_before_debatch.png"
ggsave(filename = file_path_clustering, plot = p, width = 15, height = 10)
cat("Hierarchical clustering plot saved to:", file_path_clustering, "\n")

library(sva)
exp_all_combat <- ComBat(exp_all, batch = group_list$dataset)

output_path_debatch <- "C:/R/C_Project/debatch_data/debatch_mRNA_50%union.csv"
write.csv(exp_all_combat, file = output_path_debatch, row.names = TRUE)
cat("Debatch data saved to:", output_path_debatch, "\n")

group_list <- data.frame(
  sample  = colnames(exp_all_combat),
  dataset = c(
    rep("GSE30784",   ncol(exp1)),
    rep("GSE41613",   ncol(exp2)),
    rep("GSE65858",   ncol(exp3)),
    rep("GSE75538",   ncol(exp4)),
    rep("GSE85446",   ncol(exp5)),
    rep("TCGA_HNSC",  ncol(exp6))
  )
)
rownames(group_list) <- group_list$sample

library(umap)
umap_result_combat <- umap(t(exp_all_combat))
umap_data_combat <- data.frame(
  UMAP1   = umap_result_combat$layout[, 1],
  UMAP2   = umap_result_combat$layout[, 2],
  dataset = group_list$dataset
)

umap_plot_combat <- ggplot(umap_data_combat, aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Gene Expression Data After Batch Effect Removal",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

output_path_combat <- "C:/R/C_Project/debatch_data/UMAP_plot_after_debatch.png"
ggsave(output_path_combat, plot = umap_plot_combat, width = 8, height = 6)
cat("UMAP plot after batch effect removal saved to:", output_path_combat, "\n")

pca_result <- prcomp(t(exp_all_combat), scale. = TRUE)
pca_data   <- data.frame(
  PC1     = pca_result$x[, 1],
  PC2     = pca_result$x[, 2],
  dataset = group_list$dataset
)
library(ggplot2)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "PCA of Gene Expression Data", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()
output_path <- "C:/R/C_Project/debatch_data/PCA_plot_after_debatch.png"
ggsave(output_path, plot = pca_plot, width = 8, height = 6)
cat("PCA plot saved to:", output_path, "\n")

library(tidyverse)
exp_all_long <- as.data.frame(exp_all_combat) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "expression", -gene)
exp_all_long$batch <- case_when(
  grepl("^GSE85446", exp_all_long$sample) ~ "GSE85446",
  grepl("^GSE41613", exp_all_long$sample) ~ "GSE41613",
  grepl("^GSE65858", exp_all_long$sample) ~ "GSE65858",
  grepl("^GSE75538", exp_all_long$sample) ~ "GSE75538",
  grepl("^GSE30784", exp_all_long$sample) ~ "GSE30784",
  grepl("^TCGA_HNSC", exp_all_long$sample) ~ "TCGA_HNSC"
)
dataset_colors <- c(
  "GSE30784" = rgb(248, 118, 109, maxColorValue = 255),
  "GSE41613" = rgb(183, 159, 0, maxColorValue = 255),
  "GSE65858" = rgb(0, 186, 56, maxColorValue = 255),
  "GSE75538" = rgb(0, 191, 196, maxColorValue = 255),
  "GSE85446" = rgb(97, 156, 255, maxColorValue = 255),
  "TCGA_HNSC" = rgb(245, 100, 227, maxColorValue = 255)
)
exp_all_long$batch <- factor(exp_all_long$batch, levels = names(dataset_colors))
p <- ggplot(exp_all_long, aes(x = sample, y = expression, fill = batch, colour = batch)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  scale_colour_manual(values = dataset_colors) +
  theme_bw() +
  labs(title = "Boxplot of Gene Expression by Sample",
       x = "Sample",
       y = "Gene Expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
file_path <- "C:/R/C_Project/debatch_data/boxplot_after_debatch.png"
ggsave(filename = file_path, plot = p, width = 15, height = 10)
cat("Boxplot saved to:", file_path, "\n")

library(tidyverse)
library(ggdendro)
data_to_cluster <- exp_all_combat
tmp_data <- dist(as.data.frame(t(data_to_cluster)))
clustering <- hclust(tmp_data, method = "average")
dendro_data <- dendro_data(clustering)
dendro_data$labels$batch <- factor(group_list$dataset[match(dendro_data$labels$label, group_list$sample)])
dataset_colors <- c(
  "GSE30784" = rgb(248, 118, 109, maxColorValue = 255),
  "GSE41613" = rgb(183, 159, 0, maxColorValue = 255),
  "GSE65858" = rgb(0, 186, 56, maxColorValue = 255),
  "GSE75538" = rgb(0, 191, 196, maxColorValue = 255),
  "GSE85446" = rgb(97, 156, 255, maxColorValue = 255),
  "TCGA_HNSC" = rgb(245, 100, 227, maxColorValue = 255)
)
exp_all_long$batch <- factor(exp_all_long$batch, levels = names(dataset_colors))

p <- ggplot(exp_all_long, aes(x = sample, y = expression, fill = batch, colour = batch)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  scale_colour_manual(values = dataset_colors) +
  theme_bw() +
  labs(title = "Boxplot of Gene Expression by Sample",
       x = "Sample",
       y = "Gene Expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

file_path <- "C:/R/C_Project/debatch_data/boxplot_after_debatch.png"
ggsave(filename = file_path, plot = p, width = 15, height = 10)
cat("Boxplot saved to:", file_path, "\n")

library(tidyverse)
library(ggdendro)

data_to_cluster <- exp_all_combat
tmp_data <- dist(as.data.frame(t(data_to_cluster)))
clustering <- hclust(tmp_data, method = "average")
dendro_data <- dendro_data(clustering)
dendro_data$labels$batch <- factor(group_list$dataset[match(dendro_data$labels$label, group_list$sample)])

dataset_colors <- c(
  "GSE30784" = rgb(248, 118, 109, maxColorValue = 255),
  "GSE41613" = rgb(183, 159, 0, maxColorValue = 255),
  "GSE65858" = rgb(0, 186, 56, maxColorValue = 255),
  "GSE75538" = rgb(0, 191, 196, maxColorValue = 255),
  "GSE85446" = rgb(97, 156, 255, maxColorValue = 255),
  "TCGA_HNSC" = rgb(245, 100, 227, maxColorValue = 255)
)

p <- ggplot() +
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dendro_data$labels, aes(x = x, y = y, color = batch), size = 3) +
  scale_color_manual(values = dataset_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Hierarchical Clustering of Samples") +
  coord_cartesian(ylim = c(0, max(dendro_data$segments$y) * 0.6))

file_path_clustering <- "C:/R/C_Project/debatch_data/hierarchical_clustering_after_debatch.png"
ggsave(filename = file_path_clustering, plot = p, width = 15, height = 10)
cat("Hierarchical clustering plot saved to:", file_path_clustering, "\n")
