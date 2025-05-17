install_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage"
library(Seurat, lib.loc = install_path)
packageVersion("Seurat")
library(data.table)
library(tidyverse)
library(harmony)
library(dplyr)
library(clustree)
library(Matrix)
library(irlba)

input_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/Merged_NonsoupX_FilteredCounts.rds"
raw_counts <- readRDS(input_file)
print(dim(raw_counts))
print(class(raw_counts))

HNSC_MERGE <- CreateSeuratObject(
  counts = raw_counts,
  project = "Merged_NonsoupX_FilteredCounts"
)

batch_info_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/batch_information.TSV"
batch_info <- read.table(batch_info_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
batch_info$Column_Name <- as.character(batch_info$Column_Name)
HNSC_MERGE@meta.data$batch <- batch_info$Dataset_Name[
  match(rownames(HNSC_MERGE@meta.data), batch_info$Column_Name)
]
HNSC_MERGE@meta.data$group <- HNSC_MERGE@meta.data$orig.ident

HNSC_MERGE <- NormalizeData(HNSC_MERGE)
HNSC_MERGE <- FindVariableFeatures(HNSC_MERGE, selection.method = "vst", nfeatures = 5000)
HNSC_MERGE <- ScaleData(HNSC_MERGE, features = rownames(HNSC_MERGE))
HNSC_MERGE <- RunPCA(HNSC_MERGE, features = VariableFeatures(object = HNSC_MERGE), reduction.name = "pca")
HNSC_MERGE <- RunHarmony(
  HNSC_MERGE,
  reduction = "pca",
  group.by.vars = "batch",
  reduction.save = "harmony",
  plot_convergence = FALSE,
  max_iter = 200,
  theta = 3,
  dims = 1:50
)

saveRDS(
  HNSC_MERGE,
  file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/HarmonyResults_OldSeuratforAnndata.rds"
)

library(MuDataSeurat)
library(fs)
library(rhdf5)
library(SeuratDisk)

seurat_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/HarmonyResults_OldSeuratforAnndata.rds"
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1"

if (!dir_exists(output_dir)) {
  dir_create(output_dir)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory already exists:", output_dir, "\n")
}

seu <- readRDS(seurat_path)
seu <- DietSeurat(
  seu,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = rownames(seu),
  assays = "RNA",
  dimreducs = c("pca", "umap", "harmony"),
  graphs = c("RNA_nn", "RNA_snn"),
  misc = TRUE
)

i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

base_filename <- "Merged_NonsoupX_Filtered_Harmony_Counts"
h5ad_output <- file.path(output_dir, paste0(base_filename, ".h5ad"))
MuDataSeurat::WriteH5AD(seu, h5ad_output, assay = "RNA")
h5mu_output <- file.path(output_dir, paste0(base_filename, ".h5mu"))
MuDataSeurat::WriteH5MU(seu, h5mu_output)
cat("Files saved to:", output_dir, "\n")
os.chdir('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1')
file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/Merged_NonsoupX_Filtered_Harmony_Counts.h5ad"
save_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1"
adata = sc.read(file_path)
adata.raw = adata
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print(adata.obs[['pct_counts_mt']].head())
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='highly_variable_genes.png', show=False)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
ho = hm.run_harmony(adata.obsm['X_harmony'], adata.obs, 'batch', theta=3.00, max_iter_harmony=200)
adata.obsm['X_harmony'] = ho.Z_corr.T
ho = hm.run_harmony(adata.obsm['X_harmony'], adata.obs, 'group', theta=3.00, max_iter_harmony=200)
adata.obsm['X_harmony'] = ho.Z_corr.T
n_neighbors = 10
resolution = 0.4
n_pcs = 20
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_harmony')
leiden_key = f'leiden_n{n_neighbors}_r{resolution:.1f}'
sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key)
sc.tl.umap(adata)
umap_plot_filename = f'umap_leiden_n{n_neighbors}_r{resolution:.1f}_pcs{n_pcs}.png'
sc.pl.umap(
    adata, 
    color=[leiden_key], 
    title=f'UMAP Leiden n{n_neighbors} r{resolution:.1f}', 
    save=umap_plot_filename,
    show=False,  
    legend_loc='on data', 
    legend_fontsize=10,  
    legend_fontweight='bold', 
    size=8,  
    alpha=0.4
)
print(f"Completed clustering and visualization for resolution={resolution:.1f}")
final_h5ad_filename = 'MNFHC_Leiden_p20n10r0.4.h5ad'
adata.write(final_h5ad_filename)
print(f"Saved all results to {final_h5ad_filename}")
    
import scanpy as sc
import os
import shutil

working_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1"
os.chdir(working_dir)
input_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_p20n10r0.4.h5ad"
output_dir = os.path.join(working_dir, "SF1B")
figure_dir = os.path.join(output_dir, "Figure")
output_file = os.path.join(output_dir, "MNFHC_Leiden_DiffGenes.h5ad")

os.makedirs(output_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

adata = sc.read_h5ad(input_file)

sc.tl.rank_genes_groups(adata, 'leiden_n10_r0.4', method='wilcoxon')

output_dir = os.path.dirname(output_file)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

adata.write(output_file)

print(f"差异基因分析完成，结果保存为: {output_file}")
print(f"可视化图表已保存到目录: {figure_dir}")
import scanpy as sc
import numpy as np
import pandas as pd

results_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1B/MNFHC_Leiden_DiffGenes.h5ad"
output_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1B/MNFHC_Leiden_Marker.h5ad"

num = 2

adata = sc.read(results_file)

print("Available columns in adata.obs:", adata.obs.columns)

if 'rank_genes_groups' not in adata.uns:
    raise ValueError("差异基因分析结果未找到，请检查输入文件。")

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
selected_genes = []

for group in groups:
    genes = result['names'][group]
    selected_genes.extend(genes[:num])

selected_genes = list(set(selected_genes))
selected_genes = [gene for gene in selected_genes if not (gene.startswith('MT') or gene.startswith('RP'))]

if adata.raw is not None:
    valid_genes = [gene for gene in selected_genes if gene in adata.raw.var_names]
else:
    valid_genes = [gene for gene in selected_genes if gene in adata.var_names]

if not valid_genes:
    raise ValueError("未找到有效的标志基因，请检查基因名称是否匹配数据。")

print(f"有效的 marker 基因数量: {len(valid_genes)}")

adata.uns['marker_genes'] = valid_genes

result = adata.uns['rank_genes_groups']
diff_genes_df = pd.DataFrame({
    group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']
})

print("前 6 行 6 列的差异基因结果表:")
print(diff_genes_df.iloc[0:6, 0:6])

adata.write(output_file)

print(f"差异基因分析完成，结果保存为: {output_file}")
import os
work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1"
os.chdir(work_dir)
print("当前工作路径:", os.getcwd())

import scanpy as sc
input_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1B/MNFHC_Leiden_Marker.h5ad"
adata = sc.read(input_file)

if 'marker_genes' in adata.uns:
    marker_genes = adata.uns['marker_genes']
    print(f"Found marker_genes: {marker_genes}")
else:
    raise ValueError("marker_genes not found in the AnnData object.")

import matplotlib.pyplot as plt
with plt.rc_context({'font.size': 16}):
    sc.pl.dotplot(
        adata,
        marker_genes,
        groupby='leiden_n10_r0.4',
        dendrogram=True,
        standard_scale='var',
        save='.dotplot.png'
    )

dotplot_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/figures/dotplot_.dotplot.png"
new_dotplot_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/figures/SF1B.png"

try:
    os.rename(dotplot_file, new_dotplot_file)
    print(f"文件已重命名: {dotplot_file} -> {new_dotplot_file}")
except FileNotFoundError:
    print(f"文件未找到: {dotplot_file}")
