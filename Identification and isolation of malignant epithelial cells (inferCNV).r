target_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD"

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

setwd(target_dir)
getwd()

install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")

library(phylogram)
library(gridExtra)
library(grid)
library(dendextend)
library(ggthemes)
library(tidyverse)
library(gplots)
library(ggplot2)
library(Matrix)
library(dplyr)
library(infercnv)
library(miscTools)

seurat_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/HarmonyResults_OldSeuratforAnndata.rds"
leiden_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1A/leiden_clusters.tsv"
output_file <- file.path(target_dir, "MNFHC_Leiden.rds")

if (!file.exists(seurat_file)) {
  stop("Seurat 文件不存在，请检查路径！")
}

if (!file.exists(leiden_file)) {
  stop("Leiden 聚类结果文件不存在，请检查路径！")
}

seurat_obj <- tryCatch(
  readRDS(seurat_file),
  error = function(e) stop("加载 Seurat 对象时出错：", e$message)
)

library(readr)
leiden_data <- read_tsv(leiden_file, col_names = c("cell", "leiden"))

if (ncol(leiden_data) < 2) {
  stop("Leiden 聚类结果文件应包含至少两列：细胞名称和聚类标签！")
}

leiden_data$cell <- gsub("LN_", "", leiden_data$cell)
leiden_data$cell <- gsub("Dis_", "", leiden_data$cell)

missing_cells <- setdiff(Cells(seurat_obj), leiden_data$cell)
if (length(missing_cells) > 0) {
  warning(paste(
    length(missing_cells),
    "个 Seurat 对象中的细胞在 Leiden 聚类结果中找不到对应的聚类标签。"
  ))
}

seurat_obj$leiden <- leiden_data$leiden[match(Cells(seurat_obj), leiden_data$cell)]

leiden_to_cell_type <- c(
  "T Cell 01", "Macrophage", "Natural Killer Cell", "Epithelial Cell 01",
  "Monocyte", "Plasma Cell 01", "Epithelial Cell 02",
  "Cancer-associated Fibroblast 01", "Epithelial Cell 03", "B cell",
  "Endothelial Cell", "T Cell 02", "Cancer-associated Fibroblast 02",
  "Epithelial Cell 04", "T Cell 03", "Cancer-associated Fibroblast 03",
  "Plasmacytoid Dendritic Cell", "Epithelial Cell 05", "Salivary cell",
  "Mast Cell", "Langerhans Cell", "CD8-positive, Alpha-Beta T Cell",
  "Plasma Cell 02", "Smooth Muscle Cell", "Mononuclear Macrophage",
  "Ciliated Cell"
)

seurat_obj$leiden <- as.numeric(seurat_obj$leiden)
seurat_obj$cell_type <- leiden_to_cell_type[seurat_obj$leiden + 1]

cat("前几行 Leiden 聚类标签：\n"); print(head(seurat_obj$leiden))
cat("前几行 对应的 cell_type：\n"); print(head(seurat_obj$cell_type))

seurat_obj <- DietSeurat(seurat_obj, assays = "RNA", scale.data = FALSE)
saveRDS(seurat_obj, file = output_file)
cat("Seurat 对象已保存至：", output_file, "\n")

seurat_obj <- readRDS(file.path(target_dir, "MNFHC_Leiden.rds"))
head(seurat_obj@meta.data)

condition_cells <- which(
  seurat_obj@meta.data$leiden == 0 |
    grepl("Epithelial", seurat_obj@meta.data$cell_type, ignore.case = TRUE)
)

num_condition_cells <- length(condition_cells)
num_total_cells <- ncol(seurat_obj)
proportion <- num_condition_cells / num_total_cells

cat("符合条件的细胞数目：", num_condition_cells, "\n")
cat("占总细胞数的比例：", proportion, "\n")

seurat_obj_filtered <- seurat_obj[, condition_cells]
saveRDS(
  seurat_obj_filtered,
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/MNFHC_Leiden_EPC.rds"
)
cat("新的 Seurat 对象已保存为：/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/MNFHC_Leiden_EPC.rds\n")

seurat_obj <- readRDS(
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/MNFHC_Leiden_EPC.rds"
)

table(Idents(seurat_obj))
table(seurat_obj@meta.data$leiden)
table(seurat_obj@meta.data$orig.ident)

dat <- GetAssayData(seurat_obj, layer = 'counts', assay = 'RNA')
dat[1:4, 1:4]

groupinfo <- data.frame(
  v1 = colnames(dat),
  v2 = seurat_obj@meta.data$leiden
)

library(AnnoProbe)
geneInfor <- annoGene(rownames(dat), "SYMBOL", "human")
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
geneInfor <- geneInfor[with(geneInfor, order(chr_num, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]

dat <- dat[rownames(dat) %in% geneInfor[, 1], ]
dat <- dat[match(geneInfor[, 1], rownames(dat)), ]

dim(dat); head(groupinfo); dat[1:4, 1:4]
table(groupinfo$v2); dim(groupinfo)

colnames(dat) <- gsub("-", "_", colnames(dat))
expFile <- "expFile_sparse.rds"
saveRDS(dat, file = expFile)

groupFiles <- 'groupFiles.txt'
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(
  groupinfo,
  file = groupFiles,
  sep = '\t',
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

geneFile <- 'geneFile.txt'
write.table(
  geneInfor,
  file = geneFile,
  sep = '\t',
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expFile,
  annotations_file = groupFiles,
  delim = "\t",
  gene_order_file = geneFile,
  ref_group_names = c("0")
)

infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv",
  cluster_by_groups = FALSE,
  hclust_method = "ward.D2",
  analysis_mode = "samples",
  denoise = TRUE,
  HMM = FALSE,
  plot_steps = FALSE,
  leiden_resolution = "auto",
  num_threads = 20
)

rm(list = ls())
options(stringsAsFactors = FALSE)
target_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD"
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}
setwd(target_dir)
getwd()

file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/MNFHC_Leiden_EPC.rds"
seurat_obj <- readRDS(file_path)
infer_CNV_obj <- readRDS('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/infercnv/run.final.infercnv_obj')
expr <- infer_CNV_obj@expr.data
data_cnv <- as.data.frame(expr)

colnames(seurat_obj) <- gsub("-", "_", colnames(seurat_obj))
meta <- seurat_obj@meta.data

tmp1 <- expr[, infer_CNV_obj@reference_grouped_cell_indices[["0"]]]
tmp <- tmp1
down <- mean(rowMeans(tmp)) - 2 * mean(apply(tmp, 1, sd))
up <- mean(rowMeans(tmp)) + 2 * mean(apply(tmp, 1, sd))
oneCopy <- up - down

a1 <- down - 2 * oneCopy
a2 <- down - 1 * oneCopy
a3 <- up + 1 * oneCopy
a4 <- up + 2 * oneCopy

cnv_score_table <- infer_CNV_obj@expr.data
cnv_score_mat <- as.matrix(cnv_score_table)

cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A"
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B"
cnv_score_table[cnv_score_mat >= down & cnv_score_mat < up] <- "C"
cnv_score_table[cnv_score_mat >= up & cnv_score_mat <= a3] <- "D"
cnv_score_table[cnv_score_mat > a3 & cnv_score_mat <= a4] <- "E"
cnv_score_table[cnv_score_mat > a4] <- "F"

cnv_score_table_pts <- cnv_score_mat
rm(cnv_score_mat)

cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2

cell_scores_CNV <- as.data.frame(colMeans(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)

seurat_obj$cnv_score <- cell_scores_CNV[match(colnames(seurat_obj), rownames(cell_scores_CNV)), "cnv_score"]

output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(seurat_obj, file = file.path(output_dir, "MNFHC_inferCNV_EPC_with_immunecontrol.rds"))

seurat_obj_filtered <- subset(seurat_obj, subset = leiden != 0)
cat("删除 leiden 为 0 后的细胞数量：", ncol(seurat_obj_filtered), "\n")
saveRDS(seurat_obj_filtered, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/MNFHC_inferCNV_EPC.rds")
cat("Seurat 对象已保存为 MNFHC_inferCNV_EPC.rds\n")

import numpy as np
import os
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import shutil
import harmonypy as hm  

work_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E"
os.chdir(work_path)
adata = sc.read('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene_celltype.h5ad')

selected_cells = adata[adata.obs['leiden_n10_r0.4'].isin(['3', '6', '8', '17', '13'])]

adata = selected_cells

print(adata.obs)

cell_info_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/cell_info.csv'
cell_info = pd.read_csv(cell_info_path)

if 'cell_name' not in cell_info.columns:
    raise ValueError("CSV 文件中缺少 'cell_name' 列。")

cell_info['cell_name'] = cell_info['cell_name'].apply(lambda x: x.replace('-', '_'))
cell_info.set_index('cell_name', inplace=True)
cell_info['cnv_score'] = pd.to_numeric(cell_info['cnv_score'], errors='coerce')

adata.obs['temp_index'] = adata.obs.index.str.replace('-', '_')
adata.obs['temp_index'] = adata.obs['temp_index'].str.replace('Dis_', '').str.replace('LN_', '')

matching_cells = pd.Index(adata.obs['temp_index']).intersection(cell_info.index)
print(f"匹配到的细胞数量: {len(matching_cells)}")

cell_info_matched = cell_info.loc[matching_cells, ['cnv_score']]
cnv_score_mapping = cell_info_matched['cnv_score'].to_dict()

adata.obs['cnv_score'] = adata.obs['temp_index'].map(cnv_score_mapping)
adata.obs['cnv_score'] = pd.to_numeric(adata.obs['cnv_score'], errors='coerce')

adata.obs.drop(columns=['temp_index'], inplace=True)

missing_cnv_score = adata.obs['cnv_score'].isnull().sum()
total_missing = missing_cnv_score

if total_missing > 0:
    print(f"有 {total_missing} 个缺失值未被填充。")
    print(f"缺失的 cnv_score 数量: {missing_cnv_score}")
    print("请检查 cell_info.csv 和 adata.obs 的匹配情况。")
else:
    print("所有缺失值均已填充。")
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch', theta=2.0)

n_neighbors = 50
resolution = 0.2
n_pcs = 10

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_harmony')
sc.tl.leiden(
    adata,
    resolution=resolution,
    key_added='leiden_EPC_status',
    flavor='igraph',
    n_iterations=2,
    directed=False
)
sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color='leiden_EPC_status',
    title='UMAP of EPC Cells Colored by CNV Score',
    save='cnv_estimation.png',
    show=False,
    size=15,
    alpha=0.4
)

cnv_5 = adata.obs['cnv_score'].quantile(0.05)
cnv_95 = adata.obs['cnv_score'].quantile(0.95)

sc.pl.umap(
    adata,
    color='cnv_score',
    title='UMAP of EPC Cells Colored by CNV Score',
    save='umap_epc_cnv_score.png',
    show=False,
    vmin=cnv_5,
    vmax=cnv_95,
    size=15,
    alpha=0.4
)
adata.obs['EPC_status'] = pd.Categorical(
    ['None_Malignant' if x == '4' else 'Malignant' for x in adata.obs['leiden_EPC_status']],
    categories=['None_Malignant', 'Malignant']
)
adata.obs = adata.obs.drop('leiden_EPC_status', axis=1)
print(adata.obs)

import shutil

source_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/figures/umapumap_epc_cnv_score.png"
destination_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/figures/SF1E_1.png"
shutil.move(source_path, destination_path)
print(f"文件已成功转移并重命名为: {destination_path}")

file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/figures/umapcnv_estimation.png'
if os.path.exists(file_path):
    os.remove(file_path)
    print(f"文件 {file_path} 已成功删除。")
else:
    print(f"文件 {file_path} 不存在。")

output_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/MNFHC_allgene_Leiden_Marker_EPC_CNV.h5ad'
adata.write(output_path)
print(f"已将更新后的 AnnData 对象保存到 {output_path}")

import pandas as pd

columns_to_extract = ['group', 'cell_type', 'leiden_n10_r0.4', 'cnv_score', 'EPC_status']
adata_obs = adata.obs[columns_to_extract].copy()
adata_obs['cell_name'] = adata.obs.index
output_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/File4.infercnv.csv'
adata_obs.to_csv(output_file, index=False)
print(f"File saved to: {output_file}")

install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")

library(phylogram)
library(gridExtra)
library(grid)
library(dendextend)
library(ggthemes)
library(tidyverse)
library(gplots)
library(ggplot2)
library(Matrix)
library(dplyr)
library(infercnv)
library(miscTools)

dir_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/png"
if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
setwd(dir_path)
getwd()

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/MNFHC_Leiden_EPC.rds")

table(Idents(seurat_obj))
table(seurat_obj@meta.data$leiden)
table(seurat_obj@meta.data$orig.ident)

clusters <- unique(seurat_obj@meta.data$leiden)
sampled_cells <- NULL
for (cluster in clusters) {
  cluster_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$leiden == cluster, ])
  num_samples <- ceiling(length(cluster_cells) / 5)
  sampled_cluster_cells <- sample(cluster_cells, num_samples)
  sampled_cells <- c(sampled_cells, sampled_cluster_cells)
}
# 此处仅以1/5细胞进行计算，用于生成代表性图形，以避免内存过大无法生成图形的bug
seurat_obj_sampled <- subset(seurat_obj, cells = sampled_cells)
print(seurat_obj_sampled)

dat <- GetAssayData(seurat_obj_sampled, layer = 'counts', assay = 'RNA')
dat[1:4, 1:4]

groupinfo <- data.frame(v1 = colnames(dat), v2 = seurat_obj_sampled@meta.data$leiden)

library(AnnoProbe)
geneInfor <- annoGene(rownames(dat), "SYMBOL", "human")
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
geneInfor <- geneInfor[with(geneInfor, order(chr_num, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]

dat <- dat[rownames(dat) %in% geneInfor[, 1], ]
dat <- dat[match(geneInfor[, 1], rownames(dat)), ]

dim(dat)
head(groupinfo)
dat[1:4, 1:4]
table(groupinfo$v2)
dim(groupinfo)

colnames(dat) <- gsub("-", "_", colnames(dat))
expFile <- "expFile_sparse.rds"
saveRDS(dat, file = expFile)

groupFiles <- 'groupFiles.txt'
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(groupinfo, file = groupFiles, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

geneFile <- 'geneFile.txt'
write.table(geneInfor, file = geneFile, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expFile,
  annotations_file = groupFiles,
  delim = "\t",
  gene_order_file = geneFile,
  ref_group_names = c("0")
)

infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv_png",
  cluster_by_groups = FALSE,
  hclust_method = "ward.D2",
  analysis_mode = "samples",
  denoise = TRUE,
  HMM = FALSE,
  plot_steps = FALSE,
  leiden_resolution = "auto",
  num_threads = 20
)
rm(list = ls())
options(stringsAsFactors = FALSE)
source_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/png/infercnv_png/infercnv.png"
destination_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/figures/SF1C.png"
file.copy(source_path, destination_path)
cat("文件已复制并重命名为", destination_path)
folder_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1CD/png/infercnv_png"
files_to_keep <- c("run.final.infercnv_obj","infercnv.png","infercnv.preliminary.png")
all_files <- list.files(folder_path, recursive = TRUE, full.names = TRUE)
all_files_names <- basename(all_files)
files_to_delete <- all_files[!(all_files_names %in% files_to_keep)]
file.remove(files_to_delete)
cat("已删除的文件:\n")
print(files_to_delete)

