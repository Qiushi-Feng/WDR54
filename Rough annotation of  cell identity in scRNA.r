working_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1A"
os.makedirs(working_dir, exist_ok=True)

input_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1B/MNFHC_Leiden_DiffGenes.h5ad"
output_file = os.path.join(working_dir, "Differential_Genes_List.tsv")

adata = sc.read_h5ad(input_file)

result = adata.uns['rank_genes_groups']

gene_names = result['names']
gene_scores = result['scores']

diff_genes_df = pd.DataFrame(gene_names)

for group in adata.obs['leiden_n10_r0.4'].cat.categories:
    diff_genes_df[group] = gene_names[group]

diff_genes_df.to_csv(output_file, sep='\t', index=False)

print(f"差异基因列表已保存为: {output_file}")

import scanpy as sc

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_p20n10r0.4.h5ad"
output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1A/leiden_clusters.tsv"

adata = sc.read_h5ad(file_path)

cluster_column = 'leiden_n10_r0.4'

if cluster_column in adata.obs.columns:
    leiden_clusters = adata.obs[[cluster_column]]
    leiden_clusters.to_csv(output_path, sep='\t')
    print(f"Leiden 聚类结果已保存到 {output_path}")
else:
    print(f"\n未找到名为 '{cluster_column}' 的聚类列，请检查列名。")

## 手动注释
# 注释方法
# ACT（手动注释，100个差异基因，不分先后。突出差异基因全局景观），选定口腔和颈部淋巴结，输入每个文件的前50个标志基因
# PanglaoDB（手动注释，最差异的5个基因。凸显差异基因的权重），下载所有细胞标志物，根据经典程度检索每个细胞群前5个可以检索到的基因做出判断，最终保存至少两个基因能够明确确定的细胞（如果没有，则按照经典程度保存前三个）,总之最终结果确定为3个
# Enrichr（手动注释，100个差异基因，不分先后。突出多数据库结果）, 输入前100个基因，选择前5个库中的p最显著合理解释

## 自动注释
install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")

library(readr)
library(SingleR)
library(celldex)
library(pheatmap)
library(dplyr)
library(textshape)
library(scater)
library(SingleCellExperiment)

seurat_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/HarmonyResults_OldSeuratforAnndata.rds"
leiden_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1A/leiden_clusters.tsv"
output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1A/singleR_Cellannotation_ReferGSE188737.TSV"

if (!file.exists(seurat_file)) {
  stop("Seurat 文件不存在，请检查路径！")
}
if (!file.exists(leiden_file)) {
  stop("Leiden 聚类结果文件不存在，请检查路径！")
}

seurat_obj <- tryCatch({
  readRDS(seurat_file)
}, error = function(e) {
  stop("加载 Seurat 对象时出错：", e$message)
})

leiden_data <- read_tsv(leiden_file, col_names = c("cell", "leiden"))
if (ncol(leiden_data) < 2) {
  stop("Leiden 聚类结果文件应包含至少两列：细胞名称和聚类标签！")
}

leiden_data$cell <- gsub("^LN_|^Dis_", "", leiden_data$cell)
print(head(leiden_data$cell))

missing_cells <- setdiff(Cells(seurat_obj), leiden_data$cell)
if (length(missing_cells) > 0) {
  warning(paste(length(missing_cells), "个 Seurat 对象中的细胞在 Leiden 聚类结果中找不到对应的聚类标签。"))
}

seurat_obj$leiden <- leiden_data$leiden[match(Cells(seurat_obj), leiden_data$cell)]
cat("前几行 Leiden 聚类标签：\n")
print(head(seurat_obj$leiden))

ref_obj_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Rawdata/Employed/GSE188737/GSE188737_hnscc_seu.RData"
if (!file.exists(ref_obj_path)) {
  stop("参考 Seurat 对象文件不存在，请检查路径！")
}
load(ref_obj_path)
if (!exists("hnscc_seu")) {
  stop("参考 Seurat 对象未正确加载。请检查 .RData 文件中的对象名称。")
}

hnscc_seu <- tryCatch({
  UpdateSeuratObject(hnscc_seu)
}, error = function(e) {
  stop("更新 Seurat 对象时出错：", e$message)
})
myref <- hnscc_seu

if (!"cell_type" %in% colnames(myref@meta.data)) {
  stop("'cell_type' 列不存在于参考 Seurat 对象的 meta.data 中。")
}
myref$celltype <- myref@meta.data$cell_type
cat("参考数据集中每种细胞类型的细胞数量：\n")
print(table(myref$celltype))

Idents(myref) <- "celltype"
Refassay <- log1p(AverageExpression(myref, verbose = FALSE)$RNA)
cat("Refassay 的前几行：\n")
print(head(Refassay))

ref_sce <- SingleCellExperiment(assays = list(counts = Refassay))
ref_sce <- logNormCounts(ref_sce)
cat("归一化后的数据（前4行4列）：\n")
print(logcounts(ref_sce)[1:4, 1:4])

colData(ref_sce)$Type <- colnames(Refassay)

testdata <- GetAssayData(seurat_obj, slot = "data")
pred <- SingleR(test = testdata, ref = ref_sce, labels = ref_sce$Type)
cat("SingleR 预测的细胞类型频数：\n")
print(table(pred$labels))
cat("SingleR 预测结果的前几行：\n")
print(head(pred))
cellType <- data.frame(seurat = seurat_obj$leiden, predict = pred$labels, stringsAsFactors = FALSE)
cat("Seurat 聚类结果的分布：\n")
print(sort(table(cellType$seurat)))
cat("Seurat 聚类与 SingleR 预测细胞类型的联合分布：\n")
print(table(cellType$seurat, cellType$predict))
finalmap <- cellType %>%
  group_by(seurat, predict) %>%
  summarise(Freq = n()) %>%
  arrange(seurat, desc(Freq)) %>%
  slice_head(n = 1) %>%
  arrange(seurat)
finalmap_df <- finalmap %>%
  select(seurat, predict)
cat("最终的细胞类型映射：\n")
print(finalmap_df)
write.table(finalmap_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(paste("细胞类型映射已保存为：", output_file, "\n"))
import anndata as ad

subset_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_p20n10r0.4.h5ad"
main_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/Merged_NonsoupX_Filtered_Harmony_Counts.h5ad"
output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene.h5ad"

subset = ad.read_h5ad(subset_path)
main = ad.read_h5ad(main_path)

if not subset.obs_names.equals(main.obs_names):
    raise ValueError("子集和主文件的obs名称不一致，请检查数据对齐情况。")

main.obs = subset.obs
main.obsm = subset.obsm
main.uns = subset.uns

main.write_h5ad(output_path)
print(f"更新后的文件已保存为: {output_path}")

import os
import scanpy as sc

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1"
os.chdir(work_dir)
print("当前工作路径:", os.getcwd())

input_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene.h5ad"
adata = sc.read(input_file)
adata.obs = adata.obs.drop('is_outlier', axis=1)

new_cluster_names = [
    'T Cell 01', 'Macrophage', 'Natural Killer Cell', 'Epithelial Cell 01', 'Monocyte',
    'Plasma Cell 01', 'Epithelial Cell 02', 'Cancer-associated Fibroblast 01',
    'Epithelial Cell 03', 'B cell', 'Endothelial Cell', 'T Cell 02',
    'Cancer-associated Fibroblast 02', 'Epithelial Cell 04', 'T Cell 03',
    'Cancer-associated Fibroblast 03', 'Plasmacytoid Dendritic Cell',
    'Epithelial Cell 05', 'Salivary cell', 'Mast Cell', 'Langerhans Cell',
    'CD8-positive, Alpha-Beta T Cell', 'Plasma Cell 02', 'Smooth Muscle Cell',
    'Mononuclear Macrophage', 'Ciliated Cell'
]

adata.rename_categories('leiden_n10_r0.4', new_cluster_names)
print("leiden_n10_r0.4 与细胞类型对应关系:")
print(adata.obs[['leiden_n10_r0.4']].value_counts())

sc.pl.umap(
    adata,
    color='leiden_n10_r0.4',
    legend_loc='right margin',
    legend_fontsize=10,
    legend_fontweight='bold',
    title='UMAP with Cell Types',
    frameon=False,
    save='.png',
    size=7,
    alpha=0.2
)

old_file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/umap.png"
new_file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/F1A.png"
os.rename(old_file_path, new_file_path)
print(f"文件已重命名为: {new_file_path}")

for group in adata.obs['group'].unique():
    sc.pl.umap(
        adata[adata.obs['group'] == group],
        color='leiden_n10_r0.4',
        title=f'UMAP for group {group}',
        frameon=False,
        legend_loc=None,
        size=10,
        alpha=0.4,
        save=f'F1A_{group}_umap.png'
    )

figures_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures"

for filename in os.listdir(figures_dir):
    if 'umap' in filename.lower():
        new_filename = filename.lower().replace('umap', '')
        old_file_path = os.path.join(figures_dir, filename)
        new_file_path = os.path.join(figures_dir, new_filename)
        os.rename(old_file_path, new_file_path)
        print(f"Renamed: {filename} -> {new_filename}")
