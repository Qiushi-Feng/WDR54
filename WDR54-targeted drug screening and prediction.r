## oncopredict
rm(list = ls())
library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)

query = GDCquery(project = "TCGA-HNSC", 
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "STAR - Counts")

GDCdownload(query)
dat = GDCprepare(query)
hnscc_count = assay(dat)
hnscc_count[1:4, 1:4]
save(hnscc_count, file = "C:/R/C_Project/F7/oncopredict/HNSC_count.Rdata")

rm(list = ls())
library(parallel)
cl <- makeCluster(0.75 * detectCores())
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("C:/R/C_Project/F7/oncopredict/gencode.v36.annotation.gtf.gz", format = "gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- parLapply(cl, exons_gene, function(x) { sum(width(reduce(x))) })
geneid_efflen <- data.frame(geneid = names(exons_gene_lens), efflen = as.numeric(exons_gene_lens))
dim(geneid_efflen)
save(geneid_efflen, file = "C:/R/C_Project/F7/oncopredict/gene_efflen.Rdata")
stopCluster(cl)
rm(list = ls())

library(parallel)
load("C:/R/C_Project/F7/oncopredict/gene_efflen.Rdata")
load("C:/R/C_Project/F7/oncopredict/HNSC_count.Rdata")

if (!identical(geneid_efflen$geneid, rownames(hnscc_count))) {
    stop("基因ID不匹配，请检查数据。")
}

counts2TPM <- function(count, efflength) {
    RPK <- count / (efflength / 1000)
    PMSC_rpk <- sum(RPK) / 1e6
    RPK / PMSC_rpk
}

hnscc_tpm <- apply(hnscc_count, 2, counts2TPM, efflength = geneid_efflen$efflen)
save(hnscc_tpm, file = "C:/R/C_Project/F7/oncopredict/HNSC_tpm.Rdata")
write.csv(hnscc_tpm, file = "C:/R/C_Project/F7/oncopredict/TCGA_HNSC_TPM.csv", row.names = TRUE)

library(data.table)
probe_map_file <- "C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap.csv"
tpm_file <- "C:/R/C_Project/F7/oncopredict/TCGA_HNSC_TPM.csv"

library(data.table)

probe_map_data <- fread(probe_map_file)
tpm_data <- fread(tpm_file)

probe_map_data$gene_id <- sub("\\..*", "", probe_map_data$id)
tpm_data$gene_id <- sub("\\..*", "", tpm_data$V1)

gene_map <- setNames(probe_map_data$gene, probe_map_data$gene_id)

tpm_data$V1 <- gene_map[tpm_data$gene_id]

tpm_data <- tpm_data[!is.na(tpm_data$V1) & tpm_data$V1 != "", ]

unmatched_ids <- tpm_data$gene_id[is.na(tpm_data$V1)]
cat("未被替换的基因ID数量:", length(unmatched_ids), "\n")

duplicate_genes <- tpm_data[, .N, by = V1][N > 1]
if (nrow(duplicate_genes) > 0) {
    cat("重复的基因名:\n")
    print(duplicate_genes$V1)
} else {
    cat("没有重复的基因名。\n")
}

numeric_columns <- sapply(tpm_data, is.numeric)
tpm_data_merged <- tpm_data[, lapply(.SD, mean, na.rm = TRUE), by = V1, .SDcols = names(tpm_data)[numeric_columns]]

fwrite(tpm_data_merged, tpm_file)
cat("结果已保存到:", tpm_file, "\n")

testExpr <- read.csv("C:/R/C_Project/F7/oncopredict/TCGA_HNSC_TPM.csv", row.names = 1)

total_samples <- ncol(testExpr)

expressed_counts <- rowSums(testExpr > 0)

genes_100_percent <- expressed_counts == total_samples
testExpr_clean <- testExpr[genes_100_percent, ]

cat("在 100% 样本中表达的基因数:", nrow(testExpr_clean), "\n")

write.csv(testExpr_clean, "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_TPM_clean.csv", row.names = TRUE)
cat("新的矩阵已保存到 C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_TPM_clean.csv。\n")

probe_map_data <- fread("C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap.csv")

fwrite(probe_map_data, "C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap_clean.csv")
cat("探针文件已输出为CSV格式，路径为: C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap_clean.csv\n")

expr_matrix1_file <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_fpkm_log.csv"
expr_matrix2_file <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_counts_log.csv"

probe_map_data <- fread("C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap_clean.csv")
probe_map_data$gene_id <- sub("\\..*", "", probe_map_data[[1]])
gene_map <- setNames(probe_map_data[[2]], probe_map_data$gene_id)

replace_gene_ids <- function(expr_matrix_file, output_file) {
    expr_data <- fread(expr_matrix_file)
    initial_gene_count <- nrow(expr_data)
    cat("初始基因数:", initial_gene_count, "\n")
    expr_data[[1]] <- sub("\\..*", "", expr_data[[1]])
    expr_data[[1]] <- gene_map[expr_data[[1]]]
    expr_data <- expr_data[!is.na(expr_data[[1]]) & expr_data[[1]] != "", ]
    replaced_gene_count <- nrow(expr_data)
    cat("替换后剩余的基因数:", replaced_gene_count, "\n")
    duplicate_genes <- expr_data[, .N, by = expr_data[[1]]][N > 1]
    if (nrow(duplicate_genes) > 0) {
        cat("重复的基因名:\n")
        print(duplicate_genes[[1]])
    } else {
        cat("没有重复的基因名。\n")
    }
    numeric_columns <- sapply(expr_data, is.numeric)
    expr_data <- expr_data[, lapply(.SD, mean, na.rm = TRUE), by = expr_data[[1]], .SDcols = names(expr_data)[numeric_columns]]
    final_gene_count <- nrow(expr_data)
    cat("去重并合并后剩余的基因数:", final_gene_count, "\n")
    fwrite(expr_data, output_file)
    cat("清理并保存文件:", output_file, "\n")
}
probe_map_file <- "C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap.csv"
tpm_file <- "C:/R/C_Project/F7/oncopredict/TCGA_HNSC_TPM.csv"

library(data.table)
probe_map_data <- fread(probe_map_file)
tpm_data <- fread(tpm_file)
probe_map_data$gene_id <- sub("\\..*", "", probe_map_data$id)
tpm_data$gene_id <- sub("\\..*", "", tpm_data$V1)
gene_map <- setNames(probe_map_data$gene, probe_map_data$gene_id)
tpm_data$V1 <- gene_map[tpm_data$gene_id]
tpm_data <- tpm_data[!is.na(tpm_data$V1) & tpm_data$V1 != "", ]
unmatched_ids <- tpm_data$gene_id[is.na(tpm_data$V1)]
cat("未被替换的基因ID数量:", length(unmatched_ids), "\n")
duplicate_genes <- tpm_data[, .N, by = V1][N > 1]
if (nrow(duplicate_genes) > 0) {
    cat("重复的基因名:\n")
    print(duplicate_genes$V1)
} else {
    cat("没有重复的基因名。\n")
}
numeric_columns <- sapply(tpm_data, is.numeric)
tpm_data_merged <- tpm_data[, lapply(.SD, mean, na.rm = TRUE), by = V1, .SDcols = names(tpm_data)[numeric_columns]]
fwrite(tpm_data_merged, tpm_file)
cat("结果已保存到:", tpm_file, "\n")

testExpr <- read.csv("C:/R/C_Project/F7/oncopredict/TCGA_HNSC_TPM.csv", row.names = 1)
total_samples <- ncol(testExpr)
expressed_counts <- rowSums(testExpr > 0)
genes_100_percent <- expressed_counts == total_samples
testExpr_clean <- testExpr[genes_100_percent, ]
cat("在 100% 样本中表达的基因数:", nrow(testExpr_clean), "\n")
write.csv(testExpr_clean, "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_TPM_clean.csv", row.names = TRUE)
cat("新的矩阵已保存到 C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_TPM_clean.csv。\n")

expr_matrix1_file <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_fpkm_log.csv"
expr_matrix2_file <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_counts_log.csv"

replace_gene_ids <- function(expr_matrix_file, output_file) {
    expr_data <- fread(expr_matrix_file)
    initial_gene_count <- nrow(expr_data)
    cat("初始基因数:", initial_gene_count, "\n")
    expr_data[[1]] <- sub("\\..*", "", expr_data[[1]])
    probe_map <- fread("C:/R/C_Project/F7/oncopredict/gencode.v22.annotation.gene.probeMap_clean.csv")
    probe_map$gene_id <- sub("\\..*", "", probe_map[[1]])
    gene_map <- setNames(probe_map[[2]], probe_map$gene_id)
    expr_data[[1]] <- gene_map[expr_data[[1]]]
    expr_data <- expr_data[!is.na(expr_data[[1]]) & expr_data[[1]] != "", ]
    replaced_gene_count <- nrow(expr_data)
    cat("替换后剩余的基因数:", replaced_gene_count, "\n")
    duplicate_genes <- expr_data[, .N, by = expr_data[[1]]][N > 1]
    if (nrow(duplicate_genes) > 0) {
        cat("重复的基因名:\n")
        print(duplicate_genes[[1]])
    } else {
        cat("没有重复的基因名。\n")
    }
    numeric_columns <- sapply(expr_data, is.numeric)
    expr_data <- expr_data[, lapply(.SD, mean, na.rm = TRUE), by = expr_data[[1]], .SDcols = names(expr_data)[numeric_columns]]
    final_gene_count <- nrow(expr_data)
    cat("去重并合并后剩余的基因数:", final_gene_count, "\n")
    fwrite(expr_data, output_file)
    cat("清理并保存文件:", output_file, "\n")
}

output_file1 <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_fpkm_log_clean.csv"
output_file2 <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_counts_log_clean.csv"
replace_gene_ids(expr_matrix1_file, output_file1)
replace_gene_ids(expr_matrix2_file, output_file2)

extract_expressed_genes <- function(expr_matrix_file) {
    expr_data <- fread(expr_matrix_file)
    numeric_columns <- sapply(expr_data, is.numeric)
    total_samples <- sum(numeric_columns)
    expressed_counts <- rowSums(expr_data[, .SD, .SDcols = numeric_columns] > 0)
    genes_100_percent <- expressed_counts == total_samples
    expr_data_clean <- expr_data[genes_100_percent, ]
    final_gene_count <- nrow(expr_data_clean)
    cat("在100%样本中表达的基因数:", final_gene_count, "\n")
    fwrite(expr_data_clean, expr_matrix_file)
    cat("已更新文件:", expr_matrix_file, "\n")
}

extract_expressed_genes(output_file1)
extract_expressed_genes(output_file2)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SparseArray")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
install.packages("oncoPredict")

library(oncoPredict)

GDSC2_Expr <- readRDS("C:/R/C_Project/F7/oncopredict/Training/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res  <- readRDS("C:/R/C_Project/F7/oncopredict/Training/Training Data/GDSC2_Res.rds")

dim(GDSC2_Expr)
dim(GDSC2_Res)
GDSC2_Expr[1:4, 1:4]
GDSC2_Res[1:4, 1:4]

CTRP2_Expr <- readRDS("C:/R/C_Project/F7/oncopredict/Training/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res  <- readRDS("C:/R/C_Project/F7/oncopredict/Training/Training Data/CTRP2_Res.rds")

dim(CTRP2_Expr)
dim(CTRP2_Res)
CTRP2_Expr[1:4, 1:4]
CTRP2_Res[1:4, 1:4]
setwd("/media/desk16/tly6105/Essential_Data/oncopredict")

library(oncoPredict)

CTRP2_Expr <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res  <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/CTRP2_Res.rds")

testExpr <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/Data/TCGA_HNSC_TPM_clean.csv")
rownames(testExpr) <- testExpr[, 1]
testExpr <- testExpr[, -1]

cat("CTRP2表达数据规模:", dim(CTRP2_Expr), "\n")
cat("CTRP2结果数据规模:", dim(CTRP2_Res), "\n")
cat("测试表达数据规模:", dim(testExpr), "\n")

phenotype_results <- calcPhenotype(
  trainingExprData = CTRP2_Expr,
  trainingPtype     = CTRP2_Res,
  testExprData      = as.matrix(testExpr),
  batchCorrect      = "none",
  powerTransformPhenotype = FALSE,
  minNumSamples           = 20,
  printOutput             = TRUE,
  removeLowVaryingGenes       = 0.2,
  removeLowVaringGenesFrom    = "homogenizeData"
)

res <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/calcPhenotype_Output/DrugPredictions.csv")
write.csv(res, "/media/desk16/tly6105/Essential_Data/oncopredict/Result/CTRP2_TCGA_TPM_nonlog.csv", row.names = FALSE)

dim(res)
head(res[, 1:4])


# GDSC1
setwd("/media/desk16/tly6105/Essential_Data/oncopredict")

library(oncoPredict)

GDSC1_Expr <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
GDSC1_Res  <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/GDSC1_Res.rds")

testExpr <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/Data/TCGA_HNSC_fpkm_log_clean.csv")
rownames(testExpr) <- testExpr[, 1]
testExpr <- testExpr[, -1]

cat("GDSC1表达数据规模:", dim(GDSC1_Expr), "\n")
cat("GDSC1结果数据规模:", dim(GDSC1_Res), "\n")
cat("测试表达数据规模:", dim(testExpr), "\n")

phenotype_results <- calcPhenotype(
  trainingExprData        = GDSC1_Expr,
  trainingPtype           = GDSC1_Res,
  testExprData            = as.matrix(testExpr),
  batchCorrect            = "none",
  powerTransformPhenotype = FALSE,
  minNumSamples           = 20,
  printOutput             = TRUE,
  removeLowVaryingGenes    = 0.2,
  removeLowVaringGenesFrom = "homogenizeData"
)

res <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/calcPhenotype_Output/DrugPredictions.csv")
write.csv(res, "/media/desk16/tly6105/Essential_Data/oncopredict/Result/GDSC1_TCGA_FPKM_nonlog.csv", row.names = FALSE)

dim(res)
head(res[, 1:4])

# GDSC2
setwd("/media/desk16/tly6105/Essential_Data/oncopredict")

library(oncoPredict)

GDSC2_Expr <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res  <- readRDS("/media/desk16/tly6105/Essential_Data/oncopredict/Training Data/GDSC2_Res.rds")

testExpr <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/Data/TCGA_HNSC_fpkm_log_clean.csv")
rownames(testExpr) <- testExpr[, 1]
testExpr <- testExpr[, -1]

cat("GDSC2表达数据规模:", dim(GDSC2_Expr), "\n")
cat("GDSC2结果数据规模:", dim(GDSC2_Res), "\n")
cat("测试表达数据规模:", dim(testExpr), "\n")

phenotype_results <- calcPhenotype(
  trainingExprData        = GDSC2_Expr,
  trainingPtype           = GDSC2_Res,
  testExprData            = as.matrix(testExpr),
  batchCorrect            = "none",
  powerTransformPhenotype = FALSE,
  minNumSamples           = 20,
  printOutput             = TRUE,
  removeLowVaryingGenes    = 0.2,
  removeLowVaringGenesFrom = "homogenizeData"
)

res <- read.csv("/media/desk16/tly6105/Essential_Data/oncopredict/calcPhenotype_Output/DrugPredictions.csv")
write.csv(res, "/media/desk16/tly6105/Essential_Data/oncopredict/Result/GDSC2_TCGA_FPKM_nonlog.csv", row.names = FALSE)

dim(res)
head(res[, 1:4])

ctrp2_file <- "C:/R/C_Project/F7/oncopredict/Result/CTRP2_TCGA_TPM_nonlog.csv"
tcga_file <- "C:/R/C_Project/F7/oncopredict/Data/TCGA_HNSC_TPM_clean.csv"
output_file <- "C:/R/C_Project/F7/oncopredict/Result/Analysis/CTRP2_TCGA_TPM_nonlog_carboplatin.csv"

ctrp2_data <- read.csv(ctrp2_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
carboplatin_data <- ctrp2_data[ , "carboplatin", drop = FALSE]
carboplatin_samples <- rownames(carboplatin_data)
carboplatin_samples_formatted <- gsub("\\.", "-", carboplatin_samples)
rownames(carboplatin_data) <- carboplatin_samples_formatted

tcga_data <- read.csv(tcga_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
tcga_samples_formatted <- gsub("\\.", "-", colnames(tcga_data))
colnames(tcga_data) <- tcga_samples_formatted

if ("WDR54" %in% rownames(tcga_data)) {
    common_samples <- intersect(rownames(carboplatin_data), tcga_samples_formatted)
    if (length(common_samples) > 0) {
        carboplatin_data_common <- carboplatin_data[common_samples, , drop = FALSE]
        wdr54_expression <- as.numeric(tcga_data["WDR54", common_samples])
        result_data <- data.frame(
            Sample = common_samples,
            Carboplatin = carboplatin_data_common[common_samples, "carboplatin"],
            WDR54_Expression = wdr54_expression
        )
        median_value <- median(result_data$WDR54_Expression, na.rm = TRUE)
        result_data$Group <- ifelse(result_data$WDR54_Expression <= median_value, "Low", "High")
        write.csv(result_data, file = output_file, row.names = FALSE)
        cat("分组已完成，结果已保存至：", output_file, "\n")
    } else {
        cat("没有共同的样本名，无法提取数据。\n")
    }
} else {
    cat("基因 'WDR54' 在 TCGA 数据中未找到。\n")
}
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(dplyr)
library(ggExtra)

input_file <- "C:/R/C_Project/F7/oncopredict/Result/Analysis/CTRP2_TCGA_TPM_nonlog_carboplatin.csv"
output_dir <- "C:/R/C_Project/F7/oncopredict/Result/Analysis"

data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

detect_outliers <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    return(x < lower_bound | x > upper_bound)
}

remove_outliers <- function(df, variable, group_var) {
    df_new <- data.frame()
    groups <- unique(df[[group_var]])
    for (g in groups) {
        df_group <- df[df[[group_var]] == g, ]
        outliers <- detect_outliers(df_group[[variable]])
        df_group <- df_group[!outliers, ]
        df_new <- rbind(df_new, df_group)
    }
    return(df_new)
}

clean_data <- remove_outliers(data, "Carboplatin", "Group")

clean_data$Group <- factor(clean_data$Group, levels = c("Low", "High"))

plot_half_violin_box <- function(df, variable, group_var, output_dir) {
    plot_data <- df[, c(variable, group_var)]
    colnames(plot_data) <- c("IC50", "Group")
    
    p <- ggplot(plot_data, aes(x = Group, y = IC50, fill = Group, color = Group)) +
        geom_half_violin(data = subset(plot_data, Group == "Low"), side = "l",
                         alpha = 0.2, width = 1.5, position = position_nudge(x = 0.4)) +
        geom_half_violin(data = subset(plot_data, Group == "High"), side = "r",
                         alpha = 0.2, width = 1.5, position = position_nudge(x = -0.4)) +
        stat_boxplot(geom = "errorbar", width = 0.1, color = "#4F4F4F",
                     position = position_nudge(x = c(0.35, -0.35))) +
        geom_boxplot(width = 0.1, notch = TRUE, alpha = 0.9, outlier.shape = NA,
                     position = position_nudge(x = c(0.35, -0.35))) +
        geom_quasirandom(data = subset(plot_data, Group == "Low"),
                         aes(x = as.numeric(Group) + 0.35),
                         color = "#4F4F4F", shape = 16, size = 0.5,
                         width = 0.1, alpha = 0.7) +
        geom_quasirandom(data = subset(plot_data, Group == "High"),
                         aes(x = as.numeric(Group) - 0.35),
                         color = "#4F4F4F", shape = 16, size = 0.5,
                         width = 0.1, alpha = 0.7) +
        labs(x = group_var, y = "IC50") +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_color_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
        scale_fill_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 12, color = "black", face = "bold"),
            axis.text = element_text(size = 10, color = "black"),
            text = element_text(size = 9, color = "black"),
            panel.background = element_rect(fill = alpha("#ffffff", 0.1), color = NA)
        )
    
    output_file <- file.path(output_dir, paste0(variable, "_IC50_Half_Violin_Plot.png"))
    ggsave(filename = output_file, plot = p, width = 4, height = 4, dpi = 300)
}

plot_half_violin_box(clean_data, "Carboplatin", "Group", output_dir)

library(dplyr)

read_and_filter_data <- function(input_file, group_var, variable) {
    data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
    
    detect_outliers <- function(x) {
        Q1 <- quantile(x, 0.25, na.rm = TRUE)
        Q3 <- quantile(x, 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        lower_bound <- Q1 - 1.5 * IQR
        upper_bound <- Q3 + 1.5 * IQR
        return(x < lower_bound | x > upper_bound)
    }
    
    remove_outliers <- function(df, variable, group_var) {
        df_new <- data.frame()
        groups <- unique(df[[group_var]])
        for (g in groups) {
            df_group <- df[df[[group_var]] == g, ]
            outliers <- detect_outliers(df_group[[variable]])
            df_group <- df_group[!outliers, ]
            df_new <- rbind(df_new, df_group)
        }
        return(df_new)
    }
    
    clean_data <- remove_outliers(data, variable, group_var)
    return(clean_data)
}
perform_stat_tests <- function(df, variable, group_var, output_file) {
    group_high <- df[df[[group_var]] == "High", variable]
    group_low <- df[df[[group_var]] == "Low", variable]
    t_test_result <- t.test(group_high, group_low)
    wilcox_test_result <- wilcox.test(group_high, group_low)
    analysis_results <- data.frame(
        Test = c("T-test", "Wilcoxon Rank Sum Test"),
        P_value = c(t_test_result$p.value, wilcox_test_result$p.value),
        Statistic = c(t_test_result$statistic, wilcox_test_result$statistic),
        Confidence_Interval = c(paste0(t_test_result$conf.int[1], " to ", t_test_result$conf.int[2]), NA),
        Mean_High = c(mean(group_high, na.rm = TRUE), NA),
        Mean_Low = c(mean(group_low, na.rm = TRUE), NA),
        Median_High = c(NA, median(group_high, na.rm = TRUE)),
        Median_Low = c(NA, median(group_low, na.rm = TRUE))
    )
    write.csv(analysis_results, output_file, row.names = FALSE)
}

main <- function(input_file, group_var, variable, output_dir) {
    clean_data <- read_and_filter_data(input_file, group_var, variable)
    output_file <- file.path(output_dir, paste0("CTRP2_TCGA_TPM_nonlog_", variable, "_analysis.csv"))
    perform_stat_tests(clean_data, variable, group_var, output_file)
}

input_file <- "C:/R/C_Project/F7/oncopredict/Result/Analysis/CTRP2_TCGA_TPM_nonlog_carboplatin.csv"
output_dir <- "C:/R/C_Project/F7/oncopredict/Result/Analysis"
main(input_file, "Group", "Carboplatin", output_dir)

ccle_file <- "C:/R/C_Project/F7/oncopredict/Result/GDSC1_CCLE_log.csv"
gdsc1_ccle_file <- "C:/R/C_Project/F7/oncopredict/Data/CCLE.csv"
output_dir <- "C:/R/C_Project/F7/oncopredict/Result/Analysis"
ccle_data <- read.csv(ccle_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
cisplatin_columns <- grep("^Cisplatin", colnames(ccle_data), value = TRUE)
gdsc1_ccle_data <- read.csv(gdsc1_ccle_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gdsc1_ccle_samples_formatted <- gsub("\\.", "-", colnames(gdsc1_ccle_data))
colnames(gdsc1_ccle_data) <- gdsc1_ccle_samples_formatted

if ("WDR54" %in% rownames(gdsc1_ccle_data)) {
    for (cisplatin_version in cisplatin_columns) {
        Cisplatin_data <- ccle_data[, cisplatin_version, drop = FALSE]
        Cisplatin_samples <- rownames(Cisplatin_data)
        Cisplatin_samples_formatted <- gsub("\\.", "-", Cisplatin_samples)
        rownames(Cisplatin_data) <- Cisplatin_samples_formatted
        common_samples <- intersect(rownames(Cisplatin_data), gdsc1_ccle_samples_formatted)
        if (length(common_samples) > 0) {
            Cisplatin_data_common <- Cisplatin_data[common_samples, , drop = FALSE]
            Cisplatin_data_common <- na.omit(Cisplatin_data_common)
            if (nrow(Cisplatin_data_common) > 0) {
                common_samples <- rownames(Cisplatin_data_common)
                wdr54_expression <- as.numeric(gdsc1_ccle_data["WDR54", common_samples])
                result_data <- data.frame(
                    Sample = common_samples,
                    Cisplatin = Cisplatin_data_common[common_samples, cisplatin_version],
                    WDR54_Expression = wdr54_expression
                )
                median_value <- median(result_data$WDR54_Expression, na.rm = TRUE)
                result_data$Group <- ifelse(result_data$WDR54_Expression <= median_value, "Low", "High")
                result_data$Cisplatin <- exp(result_data$Cisplatin)
                output_file <- file.path(output_dir, paste0("CCLE_GDSC1_CCLE_", cisplatin_version, ".csv"))
                write.csv(result_data, file = output_file, row.names = FALSE)
                cat("分组已完成并已将", cisplatin_version, "值逆转，结果已保存至：", output_file, "\n")
            } else {
                cat("没有有效的", cisplatin_version, "数据，所有样本值为 NA。\n")
            }
        } else {
            cat("没有共同的样本名，无法提取", cisplatin_version, "数据。\n")
        }
    }
} else {
    cat("基因 'WDR54' 在 GDSC1 CCLE 数据中未找到。\n")
}
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(dplyr)
library(ggExtra)

detect_outliers <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    return(x < lower_bound | x > upper_bound)
}

remove_outliers <- function(df, variable, group_var) {
    df_new <- data.frame()
    groups <- unique(df[[group_var]])
    for (g in groups) {
        df_group <- df[df[[group_var]] == g, ]
        outliers <- detect_outliers(df_group[[variable]])
        df_group <- df_group[!outliers, ]
        df_new <- rbind(df_new, df_group)
    }
    return(df_new)
}

plot_half_violin_box <- function(df, variable, group_var, output_file) {
    plot_data <- df[, c(variable, group_var)]
    colnames(plot_data) <- c("IC50", "Group")
    p <- ggplot(plot_data, aes(x = Group, y = IC50, fill = Group, color = Group)) +
        geom_half_violin(data = subset(plot_data, Group == "Low"), side = "l",
                         alpha = 0.2, width = 1.5, position = position_nudge(x = 0.4)) +
        geom_half_violin(data = subset(plot_data, Group == "High"), side = "r",
                         alpha = 0.2, width = 1.5, position = position_nudge(x = -0.4)) +
        stat_boxplot(geom = "errorbar", width = 0.1, color = "#4F4F4F",
                     position = position_nudge(x = c(0.35, -0.35))) +
        geom_boxplot(width = 0.1, notch = TRUE, alpha = 0.9, outlier.shape = NA,
                     position = position_nudge(x = c(0.35, -0.35))) +
        geom_quasirandom(data = subset(plot_data, Group == "Low"),
                         aes(x = as.numeric(Group) + 0.35),
                         color = "#4F4F4F", shape = 16, size = 0.5,
                         width = 0.1, alpha = 0.7) +
        geom_quasirandom(data = subset(plot_data, Group == "High"),
                         aes(x = as.numeric(Group) - 0.35),
                         color = "#4F4F4F", shape = 16, size = 0.5,
                         width = 0.1, alpha = 0.7) +
        labs(x = group_var, y = "IC50") +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_color_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
        scale_fill_manual(values = c("Low" = "#3372A6", "High" = "#AF0F11")) +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 12, color = "black", face = "bold"),
            axis.text = element_text(size = 10, color = "black"),
            text = element_text(size = 9, color = "black"),
            panel.background = element_rect(fill = alpha("#ffffff", 0.1), color = NA)
        )
    ggsave(filename = output_file, plot = p, width = 4, height = 4, dpi = 300)
}

perform_stat_tests <- function(df, variable, group_var, output_file) {
    group_high <- df[df[[group_var]] == "High", variable]
    group_low <- df[df[[group_var]] == "Low", variable]
    t_test_result <- t.test(group_high, group_low)
    wilcox_test_result <- wilcox.test(group_high, group_low)
    analysis_results <- data.frame(
        Test = c("T-test", "Wilcoxon Rank Sum Test"),
        P_value = c(t_test_result$p.value, wilcox_test_result$p.value),
        Statistic = c(t_test_result$statistic, wilcox_test_result$statistic),
        Confidence_Interval = c(paste0(t_test_result$conf.int[1], " to ", t_test_result$conf.int[2]), NA),
        Mean_High = c(mean(group_high, na.rm = TRUE), NA),
        Mean_Low = c(mean(group_low, na.rm = TRUE), NA),
        Median_High = c(NA, median(group_high, na.rm = TRUE)),
        Median_Low = c(NA, median(group_low, na.rm = TRUE))
    )
    write.csv(analysis_results, output_file, row.names = FALSE)
}

main <- function(input_dir, output_dir, file_pattern) {
    files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
    for (file in files) {
        data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        clean_data <- remove_outliers(data, "Cisplatin", "Group")
        clean_data$Group <- factor(clean_data$Group, levels = c("Low", "High"))
        file_name <- tools::file_path_sans_ext(basename(file))
        plot_output_file <- file.path(output_dir, paste0(file_name, "_IC50_Half_Violin_Plot.png"))
        plot_half_violin_box(clean_data, "Cisplatin", "Group", plot_output_file)
        analysis_output_file <- file.path(output_dir, paste0(file_name, "_analysis.csv"))
        perform_stat_tests(clean_data, "Cisplatin", "Group", analysis_output_file)
    }
}

input_dir <- "C:/R/C_Project/F7/oncopredict/Result/Analysis"
output_dir <- "C:/R/C_Project/F7/oncopredict/Result/Analysis"
file_pattern <- "CCLE_GDSC1_CCLE_Cisplatin.*\\.csv$"

main(input_dir, output_dir, file_pattern)

##cmap
input_file <- "C:/R/C_Project/F7/cmap/DEG/TCGA_HNSC_mRNA_clean.csv"
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

if ('WDR54' %in% rownames(data)) {
    WDR54_expr <- as.numeric(data['WDR54', ])
    names(WDR54_expr) <- colnames(data)
} else {
    stop("基因 'WDR54' 在数据中未找到。")
}

WDR54_expr_sorted <- sort(WDR54_expr)
n_samples <- length(WDR54_expr_sorted)
n_low <- floor(n_samples / 2)

low_samples <- names(WDR54_expr_sorted)[1:n_low]
high_samples <- names(WDR54_expr_sorted)[(n_low + 1):n_samples]

group_df <- data.frame(
    Sample = c(low_samples, high_samples),
    Group = c(rep('Low', length(low_samples)), rep('High', length(high_samples)))
)

group_df <- group_df[order(group_df$Sample), ]

output_file <- "C:/R/C_Project/F7/cmap/DEG/WDR54_grouping.csv"
write.csv(group_df, file = output_file, row.names = FALSE)

cat("分组结果已保存到：", output_file, "\n")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks")
if(!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("edgeR")) BiocManager::install("edgeR")
if(!require("limma")) BiocManager::install("limma")

if(!require("survival")) install.packages("survival")
if(!require("broom")) install.packages("broom")
if(!require("devtools")) install.packages("devtools")
if(!require("cli")) install.packages("cli")
devtools::install_github("ayueme/easyTCGA")
library(easyTCGA)

expression_file <- "C:/R/C_Project/F7/cmap/DEG/TCGA_HNSC_mRNA_clean.csv"
grouping_file <- "C:/R/C_Project/F7/cmap/DEG/WDR54_grouping.csv"
output_file <- "C:/R/C_Project/F7/cmap/DEG/TCGA_DEG.csv"

exprSet <- read.csv(expression_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
grouping <- read.csv(grouping_file, header = TRUE, stringsAsFactors = FALSE)

if (!all(grouping$Sample %in% colnames(exprSet))) {
    stop("分组文件中的样本名与表达矩阵中的样本名不一致，请检查文件。")
}

group <- as.factor(grouping$Group)
names(group) <- grouping$Sample

exprSet <- exprSet[, grouping$Sample]

diff_res <- diff_analysis(exprset = exprSet, group = group, is_count = FALSE)
diff_limma <- diff_res$deg_limma
diff_limma_filtered <- diff_limma[diff_limma$P.Value < 0.05, ]
diff_limma_sorted <- diff_limma_filtered[order(-diff_limma_filtered$logFC), ]
write.csv(diff_limma_sorted, file = output_file, row.names = FALSE)

cat("差异分析结果已保存到：", output_file, "\n")

## 将差异基因上传到cmap进行分析

library(easyTCGA)

process_expression_matrix <- function(expression_file, grouping_file, output_file, gene_of_interest = "WDR54") {
    exprSet <- read.csv(expression_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    
    if (gene_of_interest %in% rownames(exprSet)) {
        gene_expr <- as.numeric(exprSet[gene_of_interest, ])
        names(gene_expr) <- colnames(exprSet)
    } else {
        stop(paste("基因", gene_of_interest, "在数据中未找到。"))
    }
    
    gene_expr_sorted <- sort(gene_expr)
    n_samples       <- length(gene_expr_sorted)
    n_low           <- floor(n_samples / 2)
    
    low_samples     <- names(gene_expr_sorted)[1:n_low]
    high_samples    <- names(gene_expr_sorted)[(n_low + 1):n_samples]
    
    group_df <- data.frame(
        Sample = c(low_samples, high_samples),
        Group  = c(rep("Low",  length(low_samples)),
                   rep("High", length(high_samples)))
    )
    
    group_df <- group_df[order(group_df$Sample), ]
    write.csv(group_df, file = grouping_file, row.names = FALSE)
    cat("分组结果已保存到：", grouping_file, "\n")
    
    if (!all(group_df$Sample %in% colnames(exprSet))) {
        stop("分组文件中的样本名与表达矩阵中的样本名不一致，请检查文件。")
    }
    
    group  <- factor(group_df$Group, levels = c("Low", "High"))
    names(group) <- group_df$Sample
    exprSet <- exprSet[, group_df$Sample]
    
    diff_res            <- diff_analysis(exprset = exprSet, group = group, is_count = FALSE)
    diff_limma          <- diff_res$deg_limma
    diff_limma_filtered <- diff_limma[diff_limma$P.Value < 0.05, ]
    diff_limma_sorted   <- diff_limma_filtered[order(-diff_limma_filtered$logFC), ]
    
    write.csv(diff_limma_sorted, file = output_file, row.names = FALSE)
    cat("差异分析结果已保存到：", output_file, "\n")
}
base_dir <- "C:/R/C_Project/F7/cmap/DEG/"
files <- c("CDMI", "CCLE", "GSE30784_clean", "GSE41613_clean", "GSE65858_clean", "GSE75538_clean", "GSE85446_clean")

for (file_name in files) {
    expression_file <- paste0(base_dir, file_name, ".csv")
    grouping_file   <- paste0(base_dir, file_name, "_grouping.csv")
    output_file     <- paste0(base_dir, file_name, "_DEG.csv")
    process_expression_matrix(expression_file, grouping_file, output_file)
}
if (!requireNamespace("export", quietly = TRUE)) {
    install.packages("export")
}

library(data.table)
library(circlize)
library(ComplexHeatmap)
library(export)
library(dplyr)
library(tibble)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(reshape2)

input_file <- "C:/R/C_Project/F7/cmap/result/TCGA-COMPOUND.txt"
output_path <- "C:/R/C_Project/F7/cmap/result/"

res <- fread(input_file, data.table = FALSE)
res_cp <- res[res$Type == "cp", ]
res_sorted <- res_cp[order(res_cp$Score, decreasing = FALSE), ]
res_bottom50 <- res_sorted[1:50, ]
res_bottom50_selected <- res_bottom50[, c("Score", "Name", "Description")]
write.csv(res_bottom50_selected, file = paste0(output_path, "res_bottom50.csv"), row.names = FALSE)

input_raw <- res_bottom50[, c("Description", "Name")]
input <- dcast(input_raw, Description ~ Name, value.var = "Name", fun.aggregate = length)
rownames(input) <- input$Description
input <- input[, -1]
input <- as.matrix(input)

input[input > 0] <- "inhibitor"
input[input == 0] <- ""
input <- input[, order(colnames(input))]

pp <- data.frame(Name = colnames(input))
pp$Type <- "Unknown"

dd <- merge(res_bottom50[, c("Name", "Description")], pp, by = "Name", all.y = TRUE)
dd$Type[grepl('inhibitor', dd$Description, ignore.case = TRUE)] <- "Inhibitor"
dd$Type[grepl('agonist', dd$Description, ignore.case = TRUE)] <- "Agonist"
dd$Type[grepl('channel blocker', dd$Description, ignore.case = TRUE)] <- "Channel blocker"
dd$Type[grepl('exchange inhibitor', dd$Description, ignore.case = TRUE)] <- "Exchange inhibitor"
dd$Type[grepl('analog', dd$Description, ignore.case = TRUE)] <- "Analog"
dd <- dd[order(dd$Type, decreasing = FALSE), ]

alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = "white", col = "#D77071")),
  inhibitor = function(x, y, w, h) grid.points(x, y, pch = 16, size = unit(0.8, "char"))
)

ha_rowdata <- rowSums(input == "inhibitor")

type_colors <- c(
  "Agonist" = "#FDBF6F",
  "Analog" = "#33A02C",
  "Channel blocker" = "#1F78B4",
  "Exchange inhibitor" = "#6A3D9A",
  "Inhibitor" = "#23BAC5",
  "Unknown" = "#808080"
)

top_ha <- HeatmapAnnotation(
  Type = dd$Type,
  col = list(Type = type_colors),
  show_legend = TRUE,
  annotation_height = unit(30, "mm"),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 0),
  annotation_name_rot = 90
)

left_ha <- rowAnnotation(
  count = anno_barplot(-ha_rowdata, axis = FALSE, border = FALSE,
                       gp = gpar(fill = "#D77071"),
                       bar_width = 1, width = unit(2, "cm")),
  annotation_name_side = "top",
  annotation_name_rot = 0
)

column_order <- dd$Name
input <- input[, column_order, drop = FALSE]

column_names_gp <- gpar(fontsize = 10)
column_names_gp$col <- ifelse(colnames(input) == "PP-2", "red", "black")
column_names_gp$fontface <- ifelse(colnames(input) == "PP-2", "bold", "plain")

png(filename = paste0(output_path, "TCGA_oncoPrint.png"), width = 4000, height = 2000, res = 300)
oncoPrint(
  input,
  alter_fun = alter_fun,
  show_column_names = TRUE,
  column_names_side = "top",
  column_order = column_order,
  top_annotation = top_ha,
  left_annotation = left_ha,
  show_pct = FALSE,
  show_heatmap_legend = FALSE,
  column_names_gp = column_names_gp
)
decorate_annotation("Type", {
  grid.text("Mechanism of Action", unit(1, "npc") + unit(3, "mm"), just = "left")
})
dev.off()

library(data.table)
library(ggplot2)

input_file <- "C:/R/C_Project/F7/cmap/result/TCGA-COMPOUND.txt"
output_path <- "C:/R/C_Project/F7/cmap/result/"

res <- fread(input_file, data.table = FALSE)

res_cp <- res[res$Type == "cp", ]
res_cp <- res_cp[order(res_cp$Score, decreasing = TRUE), ]
res_cp$Rank <- seq_len(nrow(res_cp))

p <- ggplot(res_cp, aes(x = Rank, y = Score)) +
  geom_point(aes(color = Score), size = 3, alpha = 0.8) +
  scale_color_gradient2(low = "#CC011F", mid = "#BFBBBA", high = "#2472A3", midpoint = 0, limits = c(-100, 100)) +
  theme_minimal() +
  labs(title = "Scatter Plot of Rank vs. Score", x = "Compound Rank", y = "Score") +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")
  )

ggsave(filename = paste0(output_path, "TCGA_scatter_plot_with_ticks.png"), plot = p, width = 10, height = 7, dpi = 300)

