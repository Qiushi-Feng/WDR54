library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(SoupX)
library(DoubletFinder)
library(Matrix)
library(BiocParallel)
library(qs)

files <- list.files('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Rawcounts',
                    pattern="\\.rds$", full.names=TRUE)

for (original_file_path in files) {
  original_file_name <- basename(original_file_path)
  new_file_name <- gsub("RawCounts.rds","ProcessedCounts.rds",original_file_name)
  save_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step1_Standard_rough_processing/nonsoupx"
  new_file_path <- file.path(save_dir,new_file_name)
  log_file_path <- file.path(save_dir,gsub("\\.rds$",".log",new_file_name))
  log_con <- file(log_file_path, open="wt")
  log_message <- function(message) {
    cat(message, "\n")
    writeLines(message, log_con)
  }
  log_message(paste0("开始处理文件：", original_file_name))
  data <- readRDS(original_file_path)
  if (!inherits(data,"dgCMatrix")) {
    data <- as(data,"dgCMatrix")
    log_message("原始数据已转换为 dgCMatrix 格式。")
  }
  filtered_cells <- colSums(data>0) > 200
  filtered_genes <- rowSums(data>0) > 3
  filtered_data <- data[filtered_genes,filtered_cells]
  log_message(paste("过滤后保留的细胞数量：", sum(filtered_cells)))
  log_message(paste("过滤后保留的基因数量：", sum(filtered_genes)))
  HNSC <- CreateSeuratObject(counts=filtered_data, project="GSE227156_Primary_21240647")
  log_message("Seurat 对象已创建。")
  log_message(paste("Seurat 对象包含", dim(HNSC)[2], "个细胞和", dim(HNSC)[1], "个基因。"))
  if (!"nFeature_RNA"%in%colnames(HNSC@meta.data)) {
    HNSC$nFeature_RNA <- Matrix::colSums(GetAssayData(HNSC,assay="RNA",layer="counts")>0)
    log_message("添加 nFeature_RNA 到元数据。")
  }
  if (!"nCount_RNA"%in%colnames(HNSC@meta.data)) {
    HNSC$nCount_RNA <- Matrix::colSums(GetAssayData(HNSC,assay="RNA",layer="counts"))
    log_message("添加 nCount_RNA 到元数据。")
  }
  HNSC[['percent.mt']] <- PercentageFeatureSet(HNSC,pattern="^MT-")
  HNSC[['percent.rb']] <- PercentageFeatureSet(HNSC,pattern="^RP[SL]")
  log_message("计算了线粒体和核糖体基因的百分比。")
  if (all(c("nFeature_RNA","percent.mt","percent.rb")%in%colnames(HNSC@meta.data))) {
    HNSC <- subset(HNSC, subset=percent.mt<20 & percent.rb<50)
    log_message("细胞过滤完成。")
  } else {
    stop("元数据缺失，无法进行细胞过滤")
  }
  HNSC <- NormalizeData(HNSC, normalization.method="LogNormalize", scale.factor=10000)
  HNSC <- FindVariableFeatures(HNSC, selection.method="vst", nfeatures=5000)
  log_message("数据规范化和高度可变基因识别完成。")
  all.genes <- rownames(HNSC)
  HNSC <- ScaleData(HNSC, features=all.genes)
  HNSC <- RunPCA(HNSC, features=VariableFeatures(object=HNSC))
  log_message("数据缩放和 PCA 降维完成。")
  HNSC <- FindNeighbors(HNSC, dims=1:10)
  HNSC <- FindClusters(HNSC, resolution=0.5)
  HNSC <- RunUMAP(HNSC, dims=1:10)
  log_message("聚类和 UMAP 降维完成。")
  sweep.res.list_HNSC <- paramSweep(HNSC, PCs=1:10, sct=FALSE)
  sweep.stats_HNSC <- summarizeSweep(sweep.res.list_HNSC, GT=FALSE)
  pdf_file <- "find_pK_output.pdf"
  tryCatch({
    pdf(pdf_file, useDingbats=FALSE)
    bcmvn_HNSC <- find.pK(sweep.stats_HNSC)
    dev.off()
    if (file.exists(pdf_file)) file.remove(pdf_file)
    log_message("find.pK 可视化结果已成功保存并删除。")
  }, error=function(e) {
    if (exists("dev.off")) dev.off()
    if (file.exists(pdf_file)) file.remove(pdf_file)
    log_message(paste("find.pK 运行时出错:", e$message))
  })
  mpK <- as.numeric(as.vector(bcmvn_HNSC$pK[which.max(bcmvn_HNSC$BCmetric)]))
  log_message(paste0("最佳 pK 值为：", mpK))
  annotations <- HNSC$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  log_message(paste0("同类型双重细胞比例为：", homotypic.prop))
  nExp_poi <- round(0.075 * nrow(HNSC@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  log_message(paste0("调整后的预期双重细胞数为：", nExp_poi.adj))
  HNSC <- doubletFinder(HNSC, PCs=1:10, pN=0.25, pK=mpK, nExp=nExp_poi.adj, reuse.pANN=FALSE, sct=FALSE)
  df_classification_cols <- grep("^DF.classifications", colnames(HNSC@meta.data), value=TRUE)
  if (length(df_classification_cols)>0) {
    latest_classification_col <- df_classification_cols[length(df_classification_cols)]
    log_message(paste0("使用的双重细胞分类列名为：", latest_classification_col))
  } else {
    log_message("未找到 DoubletFinder 分类结果列。")
  }
  doublet_cells <- which(HNSC@meta.data[[latest_classification_col]]=="Doublet")
  log_message(paste0("检测到的双重细胞数：", length(doublet_cells)))
  if (length(doublet_cells)>0) {
    HNSC <- subset(HNSC, cells=setdiff(colnames(HNSC), colnames(HNSC)[doublet_cells]))
    log_message("双重细胞已删除。")
  }
  group_A <- c("PTPRC")
  group_B <- c("EPCAM","KRT19","KRT8")
  group_C <- c("LUM","DCN","COL1A1","ACTA2")
  group_D <- c("PECAM1","VWF","ENG")
  rna_assay <- HNSC@assays$RNA
  counts_matrix <- rna_assay@layers$counts
  valid_genes_A <- group_A[group_A%in%rownames(counts_matrix)]
  valid_genes_B <- group_B[group_B%in%rownames(counts_matrix)]
  valid_genes_C <- group_C[group_C%in%rownames(counts_matrix)]
  valid_genes_D <- group_D[group_D%in%rownames(counts_matrix)]
  if (length(valid_genes_A)==0) cat("警告: A 组基因未找到\n")
  if (length(valid_genes_B)==0) cat("警告: B 组基因未找到\n")
  if (length(valid_genes_C)==0) cat("警告: C 组基因未找到\n")
  if (length(valid_genes_D)==0) cat("警告: D 组基因未找到\n")
  counts_A <- counts_matrix[valid_genes_A,,drop=FALSE]
  counts_B <- counts_matrix[valid_genes_B,,drop=FALSE]
  counts_C <- counts_matrix[valid_genes_C,,drop=FALSE]
  counts_D <- counts_matrix[valid_genes_D,,drop=FALSE]
  num_groups_non_zero <- rowSums(rbind(
    colSums(counts_A!=0)>0,
    colSums(counts_B!=0)>0,
    colSums(counts_C!=0)>0,
    colSums(counts_D!=0)>0
  ))
  cells_to_remove <- colnames(counts_matrix)[num_groups_non_zero>=2]
  if (length(cells_to_remove)>0) {
    HNSC <- subset(HNSC, cells=setdiff(colnames(HNSC), cells_to_remove))
    if ("reductions"%in%slotNames(HNSC)) {
      for (reduction in names(HNSC@reductions)) {
        HNSC@reductions[[reduction]]@cell.embeddings <-
          HNSC@reductions[[reduction]]@cell.embeddings[!rownames(HNSC@reductions[[reduction]]@cell.embeddings)%in%cells_to_remove,]
      }
    }
    if ("graphs"%in%slotNames(HNSC)) {
      for (graph in names(HNSC@graphs)) {
        HNSC@graphs[[graph]] <-
          HNSC@graphs[[graph]][!rownames(HNSC@graphs[[graph]])%in%cells_to_remove,
                                !colnames(HNSC@graphs[[graph]])%in%cells_to_remove]
      }
    }
    cat("细胞已成功删除。\n")
  } else {
    cat("没有需要删除的细胞。\n")
  }
  HNSC <- NormalizeData(HNSC, normalization.method="LogNormalize", scale.factor=10000)
  HNSC <- FindVariableFeatures(HNSC, selection.method="vst", nfeatures=5000)
  all.genes <- rownames(HNSC)
  HNSC <- ScaleData(HNSC, features=all.genes)
  HNSC <- RunPCA(HNSC, features=VariableFeatures(object=HNSC))
  HNSC <- FindNeighbors(HNSC, dims=1:10)
  HNSC <- FindClusters(HNSC, resolution=0.5)
  HNSC <- RunUMAP(HNSC, dims=1:10)
  cat("数据处理完成！\n")
  cat("标准化和降维分析已完成，并覆盖原有 Seurat 对象。\n")
  saveRDS(HNSC, file=new_file_path)
  log_message(paste("已覆盖保存处理后的 Seurat 对象:", new_file_path))
  close(log_con)
  cat("日志文件已关闭。\n")
  rm(list=setdiff(ls(), c("files", "save_dir", "log_message")))
  gc()
  cat("内存已释放。\n")
}

cat("所有文件处理完成。\n")

processed_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step1_Standard_rough_processing/nonsoupx"
rawcounts_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Rawcounts"
output_file <- file.path(processed_dir, "counts_comparison_results.tsv")

processed_files <- list.files(processed_dir, pattern="\\.rds$", full.names=TRUE)

results <- data.frame(
  Processed_File=character(),
  Processed_Rows=integer(),
  Processed_Cols=integer(),
  Raw_File=character(),
  Raw_Rows=integer(),
  Raw_Cols=integer(),
  Row_Diff=integer(),
  Col_Diff=integer(),
  stringsAsFactors=FALSE
)

for (processed_file in processed_files) {
  seurat_obj <- readRDS(processed_file)
  if ("RNA"%in%names(seurat_obj@assays)) {
    rna_assay <- seurat_obj@assays$RNA
    counts_matrix <- rna_assay@layers$counts
    processed_counts_rows <- nrow(counts_matrix)
    processed_counts_cols <- ncol(counts_matrix)
  } else {
    next
  }

  pattern <- sub("_ProcessedCounts\\.rds$", "", basename(processed_file))
  raw_files <- list.files(rawcounts_dir, pattern=paste0("^", pattern), full.names=TRUE)
  if (length(raw_files)==0) next

  raw_counts_file <- raw_files[1]
  raw_counts_matrix <- readRDS(raw_counts_file)
  if (!inherits(raw_counts_matrix,"dgCMatrix")) next

  raw_counts_rows <- nrow(raw_counts_matrix)
  raw_counts_cols <- ncol(raw_counts_matrix)

  row_diff <- processed_counts_rows - raw_counts_rows
  col_diff <- processed_counts_cols - raw_counts_cols

  results <- rbind(
    results,
    data.frame(
      Processed_File=basename(processed_file),
      Processed_Rows=processed_counts_rows,
      Processed_Cols=processed_counts_cols,
      Raw_File=basename(raw_counts_file),
      Raw_Rows=raw_counts_rows,
      Raw_Cols=raw_counts_cols,
      Row_Diff=row_diff,
      Col_Diff=col_diff,
      stringsAsFactors=FALSE
    )
  )
}

write.table(results, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)

cat("分析完成，结果已保存到:", output_file, "\n")
library(Seurat)
library(Matrix)

source_directory <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step1_Standard_rough_processing/nonsoupx"
output_directory <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix/nonsoupX"

pattern <- "(GSE188737|GSE227156|GSE234933_Metastasis).*\\.rds$"
source_files <- list.files(source_directory, pattern = pattern, full.names = TRUE)
print(source_files)

for (source_file in source_files) {
  seurat_object <- readRDS(source_file)
  rna_assay <- seurat_object@assays$RNA
  cell_names <- rownames(rna_assay@cells@.Data)
  feature_names <- rownames(rna_assay@features@.Data)
  counts_matrix <- rna_assay@layers$counts
  dimnames(counts_matrix) <- list(feature_names, cell_names)
  counts_sparse <- as(counts_matrix, "CsparseMatrix")
  output_file_name <- gsub("ProcessedCounts", "FilteredCounts", basename(source_file))
  output_file_path <- file.path(output_directory, output_file_name)
  saveRDS(counts_sparse, file = output_file_path)
  cat("Sparse counts matrix saved to: ", output_file_path, "\n")
  rm(seurat_object, rna_assay, counts_matrix, counts_sparse)
  gc()
}

rds_directory <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix/nonsoupX"
if (!dir.exists(rds_directory)) stop(paste("The directory", rds_directory, "does not exist."))
rds_files <- list.files(rds_directory, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in the specified directory.")

total_columns <- 0
unique_row_names <- character(0)
shared_row_names <- NULL
wdr54_shared <- TRUE

for (file in rds_files) {
  sparse_matrix <- readRDS(file)
  if (!inherits(sparse_matrix, "dgCMatrix")) stop(paste("File", file, "is not a dgCMatrix object."))
  current_row_names <- rownames(sparse_matrix)
  total_columns <- total_columns + ncol(sparse_matrix)
  unique_row_names <- union(unique_row_names, current_row_names)
  if (is.null(shared_row_names)) {
    shared_row_names <- current_row_names
  } else {
    shared_row_names <- intersect(shared_row_names, current_row_names)
  }
  if (!"WDR54" %in% current_row_names) {
    wdr54_shared <- FALSE
  }
}

cat("Total Columns Across All Files:", total_columns, "\n")
cat("Unique Row Names Count:", length(unique_row_names), "\n")
cat("Shared Row Names Count:", length(shared_row_names), "\n")
cat("Is WDR54 Shared Across All Files:", wdr54_shared, "\n")

library(Matrix)
library(dplyr)

input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix/nonsoupX"
output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix"
rds_files <- list.files(input_path, pattern = "\\.rds$", full.names = TRUE)

list_matrices <- lapply(rds_files, function(file) {
  mat <- readRDS(file)
  return(mat)
})

gene_union <- Reduce(union, lapply(list_matrices, rownames))

list_matrices <- lapply(list_matrices, function(mat) {
  missing_genes <- setdiff(gene_union, rownames(mat))
  if (length(missing_genes) > 0) {
    missing_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(missing_matrix) <- missing_genes
    colnames(missing_matrix) <- colnames(mat)
    mat <- rbind(mat, missing_matrix)
  }
  return(mat[gene_union, , drop = FALSE])
})
num_groups <- 12
matrix_groups <- split(list_matrices, rep(1:num_groups, length.out = length(list_matrices)))

for (i in 1:num_groups) {
  combined_matrix <- do.call(cbind, matrix_groups[[i]])
  output_file <- file.path(output_path, paste0("group1_", sprintf("%02d", i), ".rds"))
  saveRDS(combined_matrix, output_file)
  cat("Saved group", i, "to", output_file, "\n")
}

rm(list = ls())
gc()

input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix"
output_path <- input_path

rds_files_group1 <- list.files(input_path, pattern = "^group1_\\d{2}\\.rds$", full.names = TRUE)

list_matrices_group1 <- lapply(rds_files_group1, function(file) {
  mat <- readRDS(file)
  return(mat)
})

gene_union_group1 <- Reduce(union, lapply(list_matrices_group1, rownames))

list_matrices_group1 <- lapply(list_matrices_group1, function(mat) {
  missing_genes <- setdiff(gene_union_group1, rownames(mat))
  if (length(missing_genes) > 0) {
    missing_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(missing_matrix) <- missing_genes
    colnames(missing_matrix) <- colnames(mat)
    mat <- rbind(mat, missing_matrix)
  }
  return(mat[gene_union_group1, , drop = FALSE])
})

num_groups2 <- 3
matrix_groups2 <- split(list_matrices_group1, rep(1:num_groups2, length.out = length(list_matrices_group1)))

for (i in 1:num_groups2) {
  combined_matrix2 <- do.call(cbind, matrix_groups2[[i]])
  output_file2 <- file.path(output_path, paste0("group2_", sprintf("%02d", i), ".rds"))
  saveRDS(combined_matrix2, output_file2)
  cat("Saved group2", i, "to", output_file2, "\n")
}

rm(list = ls())
gc()
library(Matrix)

input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix"
output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony"

rds_files_group2 <- list.files(input_path, pattern = "^group2_\\d{2}\\.rds$", full.names = TRUE)

list_matrices_group2 <- lapply(rds_files_group2, function(file) {
  mat <- readRDS(file)
  return(mat)
})

gene_union_group2 <- Reduce(union, lapply(list_matrices_group2, rownames))

list_matrices_group2 <- lapply(list_matrices_group2, function(mat) {
  missing_genes <- setdiff(gene_union_group2, rownames(mat))
  if (length(missing_genes) > 0) {
    missing_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(missing_matrix) <- missing_genes
    colnames(missing_matrix) <- colnames(mat)
    mat <- rbind(mat, missing_matrix)
  }
  return(mat[gene_union_group2, , drop = FALSE])
})

final_combined_matrix <- do.call(cbind, list_matrices_group2)

output_file_final <- file.path(output_path, "Merged_NonsoupX_FilteredCounts.rds")
saveRDS(final_combined_matrix, output_file_final)
cat("Saved final merged matrix to", output_file_final, "\n")

rm(list = ls())
gc()

target_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix"

files_to_delete <- list.files(target_path, pattern = "group.*\\.rds$", full.names = TRUE, recursive = FALSE)

if (length(files_to_delete) > 0) {
  cat("Found files to delete:\n", paste(files_to_delete, collapse = "\n"), "\n")
  file.remove(files_to_delete)
  cat("Deleted the following files:\n", paste(files_to_delete, collapse = "\n"), "\n")
} else {
  cat("No files matching the criteria were found.\n")
}

input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix/nonsoupX"
output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Merged_NonsoupX_FilteredCounts.rds"

cat("Analyzing source RDS files...\n")
rds_files <- list.files(input_path, pattern = "\\.rds$", full.names = TRUE)

source_nonzero_count <- 0
source_col_count <- 0
source_unique_rows <- c()

for (file in rds_files) {
  tryCatch({
    mat <- readRDS(file)
    if (!inherits(mat, "dgCMatrix")) stop(paste("File", file, "is not a dgCMatrix."))
    source_nonzero_count <- source_nonzero_count + length(mat@x)
    source_col_count <- source_col_count + ncol(mat)
    source_unique_rows <- union(source_unique_rows, rownames(mat))
  }, error = function(e) {
    cat("Error processing file:", file, "\n", e$message, "\n")
  })
}

source_row_count <- length(source_unique_rows)

cat("Source RDS files analyzed:\n")
cat(" - Total non-zero values:", source_nonzero_count, "\n")
cat(" - Total columns:", source_col_count, "\n")
cat(" - Unique rows:", source_row_count, "\n\n")

cat("Analyzing output file...\n")
output_matrix <- readRDS(output_file)
output_nonzero_count <- length(output_matrix@x)
output_row_count <- nrow(output_matrix)
output_col_count <- ncol(output_matrix)

cat("Output file analyzed:\n")
cat(" - Total non-zero values:", output_nonzero_count, "\n")
cat(" - Total rows:", output_row_count, "\n")
cat(" - Total columns:", output_col_count, "\n\n")

cat("Randomly selecting 100 non-zero values from output file...\n")
set.seed(42)
nonzero_indices <- which(output_matrix != 0, arr.ind = TRUE)
selected_indices <- nonzero_indices[sample(1:nrow(nonzero_indices), min(100, nrow(nonzero_indices))), ]

selected_rows <- rownames(output_matrix)[selected_indices[, 1]]
selected_cols <- colnames(output_matrix)[selected_indices[, 2]]

cat("Checking consistency across source files...\n")
discrepancy_count <- 0

for (i in seq_along(selected_rows)) {
  row_name <- selected_rows[i]
  col_name <- selected_cols[i]
  output_value <- output_matrix[row_name, col_name]
  
  found <- FALSE
  for (file in rds_files) {
    tryCatch({
      mat <- readRDS(file)
      if (row_name %in% rownames(mat) && col_name %in% colnames(mat)) {
        source_value <- mat[row_name, col_name]
        if (source_value != output_value) {
          discrepancy_count <- discrepancy_count + 1
          cat("Discrepancy found at (", row_name, ", ", col_name, "): ",
              "Source =", source_value, ", Output =", output_value, "\n")
        }
        found <- TRUE
        break
      }
    }, error = function(e) {
      cat("Error processing file:", file, "\n", e$message, "\n")
    })
  }
  if (!found) {
    cat("Value not found in source files for (", row_name, ", ", col_name, ")\n")
  }
}

cat("Consistency check completed.\n")
cat("Total discrepancies found:", discrepancy_count, "\n")
library(stringr)
input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Dispersed sparse matrix/nonsoupX"
output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/batch_information.TSV"
rds_files <- list.files(input_path, pattern="\\.rds$", full.names=TRUE)
results <- list()
for (file in rds_files) {
  file_name <- basename(file)
  dataset_name <- strsplit(file_name, "_")[[1]][1]
  data <- readRDS(file)
  column_names <- colnames(data)
  results[[file]] <- data.frame(Column_Name=column_names, Dataset_Name=dataset_name, stringsAsFactors=FALSE)
}
batch_info <- do.call(rbind, results)
write.table(batch_info, file=output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat("所有RDS文件已处理完成，结果保存为:", output_file, "\n")

library(data.table)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(harmony)
library(patchwork)
library(tidyr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(clustree)
input_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Merged_NonsoupX_FilteredCounts.rds"
figure_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/Figure/AllGeneSample"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive=TRUE)
  cat("Directory created at:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
raw_counts <- readRDS(input_file)
print(dim(raw_counts))
print(class(raw_counts))
HNSC_MERGE <- CreateSeuratObject(counts=raw_counts, project="Merged_NonsoupX_FilteredCounts")
batch_info_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/batch_information.TSV"
batch_info <- read.table(batch_info_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
batch_info$Column_Name <- as.character(batch_info$Column_Name)
HNSC_MERGE@meta.data$batch <- batch_info$Dataset_Name[match(rownames(HNSC_MERGE@meta.data), batch_info$Column_Name)]
HNSC_MERGE@meta.data$group <- HNSC_MERGE@meta.data$orig.ident
HNSC_MERGE <- NormalizeData(HNSC_MERGE)
HNSC_MERGE <- FindVariableFeatures(HNSC_MERGE, selection.method="vst", nfeatures=5000)
HNSC_MERGE <- ScaleData(HNSC_MERGE, features=rownames(HNSC_MERGE))
HNSC_MERGE <- RunPCA(HNSC_MERGE, features=VariableFeatures(object=HNSC_MERGE), reduction.name="pca")
HNSC_MERGE <- RunHarmony(
  HNSC_MERGE,
  reduction="pca",
  group.by.vars="batch",
  reduction.save="harmony",
  plot_convergence=FALSE,
  max_iter=200,
  theta=3,
  dims=1:50
)
saveRDS(HNSC_MERGE, file="/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Merged_NonsoupX_Filtered_Harmony_Counts.rds")
