import os
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import shutil
import numpy as np
import harmonypy as hm  # 导入 harmonypy

matplotlib.use('Agg')

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC"
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

os.chdir(work_dir)
print(f"Current working directory: {os.getcwd()}")

adata = sc.read('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/MNFHC_allgene_Leiden_Marker_EPC_CNV.h5ad')

category_ePC_status = adata.obs.groupby(['orig.ident', 'EPC_status']).size().unstack(fill_value=0)
category_ePC_status_percentage = category_ePC_status.div(category_ePC_status.sum(axis=1), axis=0) * 100
print(category_ePC_status_percentage)

if 'EPC_status' not in adata.obs:
    raise ValueError("EPC_status column not found in the AnnData object.")
malignant_cells = adata[adata.obs['EPC_status'] == 'Malignant', :]
print(malignant_cells.obs)

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene_celltype.h5ad"
adata = sc.read(file_path)
print(adata.obs.head())

if 'cell_type' not in adata.obs:
    raise ValueError("在 AnnData 对象中未找到 'cell_type' 列。")

total_cells = adata.shape[0]
tumor_associated_fibroblasts = adata[adata.obs['cell_type'].str.contains('Cancer-associated Fibroblast', case=False, na=False), :]
fibroblast_count = tumor_associated_fibroblasts.shape[0]
non_fibroblast_count = total_cells - fibroblast_count

print(f"总细胞数: {total_cells}")
print(f"提取的含有 'fib' 的细胞数量（tumor_associated_fibroblasts）: {fibroblast_count}")
print(f"未提取的细胞数量 (非 'fib' 细胞): {non_fibroblast_count}")
print(tumor_associated_fibroblasts.obs)

combined = malignant_cells.concatenate(
    tumor_associated_fibroblasts,
    batch_key='dataset',
    batch_categories=['MalignantEPC', 'CAF']
)
print(combined)
print(combined.obs.head())

adata = combined
original_indices = adata.obs_names.tolist()

def remove_last_dash(index):
    if '-' in index:
        return '-'.join(index.split('-')[:-1])
    else:
        return index

new_indices = [remove_last_dash(idx) for idx in original_indices]
adata.obs_names = pd.Index(new_indices)

print(adata)
print(adata.obs.head())









sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')

ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch', theta=0.5)
adata.obsm['X_harmony'] = ho.Z_corr.T

subset_cells = adata.obs['dataset'].isin(['CAF'])
adata_partial = adata[subset_cells].copy()
print(f"新对象 adata_partial 中包含 {adata_partial.n_obs} 个细胞")

ho = hm.run_harmony(adata_partial.obsm['X_harmony'], adata_partial.obs, 'batch', theta=1.0)
adata_partial.obsm['X_harmony'] = ho.Z_corr.T

cell_names_partial = adata_partial.obs_names
indices = adata.obs_names.get_indexer(cell_names_partial)
adata.obsm['X_harmony'][indices, :] = adata_partial.obsm['X_harmony']

print("成功将 adata_partial 的 X_harmony 数据写回到 adata 中对应的细胞。")

n_neighbors = 100
resolution = 0.4
n_pcs = 30

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_harmony')

leiden_mEPC_CAF = f'leiden_mEPC_CAF_n{n_neighbors}_r{resolution:.1f}'

sc.tl.leiden(
    adata,
    resolution=resolution,
    key_added=leiden_mEPC_CAF,
    flavor='igraph',
    n_iterations=2,
    directed=False
)

sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color=[leiden_mEPC_CAF],
    save=umap_plot_filename,
    show=False,
    size=20,
    alpha=0.4,
    legend_loc='on data',
    legend_fontsize=12,
    legend_fontweight='bold'
)

leiden_dataset_ratio = adata.obs.groupby(leiden_mEPC_CAF)['dataset'].value_counts(normalize=True)
print(leiden_dataset_ratio)






import pandas as pd

colors_key = f'{leiden_mEPC}_colors'

if colors_key in malignant_cells.uns:
    leiden_colors = malignant_cells.uns[colors_key]
    print(f"已提取 '{colors_key}' 配色方案。")
else:
    raise KeyError(f"'{colors_key}' 未在 malignant_cells.uns 中找到。")

cluster_proportions = pd.crosstab(
    malignant_cells.obs['orig.ident'],
    malignant_cells.obs[leiden_mEPC],
    normalize='index'
) * 100

desired_order = ["Pr", "LN_Me", "Dis_Me"]
existing_order = [ident for ident in desired_order if ident in cluster_proportions.index]
cluster_proportions_sorted = cluster_proportions.loc[existing_order]

unique_clusters_sorted = sorted(cluster_proportions_sorted.columns)

if len(leiden_colors) < len(unique_clusters_sorted):
    raise ValueError(
        f"配色方案的颜色数量少于 '{leiden_mEPC}' 聚类标签的数量。请确保配色方案包含足够的颜色。"
    )

color_mapping = {
    cluster: color for cluster, color in zip(unique_clusters_sorted, leiden_colors)
}

plt.ioff()
fig, ax = plt.subplots(figsize=(12, 8))
bottom = pd.Series([0] * len(cluster_proportions_sorted), index=cluster_proportions_sorted.index)

for cluster in unique_clusters_sorted:
    ax.barh(
        cluster_proportions_sorted.index,
        cluster_proportions_sorted[cluster],
        left=bottom,
        color=color_mapping[cluster],
        edgecolor='black',
        label=str(cluster)
    )
    bottom += cluster_proportions_sorted[cluster]

ax.set_xlabel('Proportion (%)', fontsize=14)
ax.set_ylabel('orig.ident', fontsize=14)
ax.set_title(f'Composition of {leiden_mEPC} in Different orig.ident', fontsize=16)

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles,
    labels,
    title=leiden_mEPC,
    bbox_to_anchor=(1.05, 1),
    loc='upper left'
)

ax.invert_yaxis()
plt.tight_layout()

output_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/F1D_MalignantEPC.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"图像已保存到 {output_path}")

output_file = os.path.join(work_dir, 'MNFHC_MalignantEPC.h5ad')
malignant_cells.write(output_file)
print(f"Saved AnnData object to {output_file}")
CAF_cells = adata.obs['leiden_mEPC_CAF_division'].str.contains('CAF')
tumor_associated_fibroblasts = adata[CAF_cells].copy()
print(tumor_associated_fibroblasts)

n_neighbors = 100
resolution = 0.4
n_pcs = 30
leiden_CAF = f'leidenCAF_n{n_neighbors}_r{resolution:.1f}'

sc.tl.leiden(
    tumor_associated_fibroblasts,
    resolution=resolution,
    key_added=leiden_CAF,
    flavor='igraph',
    n_iterations=2,
    directed=False
)

tumor_associated_fibroblasts.obs[leiden_CAF] = (
    'CAF ' +
    (tumor_associated_fibroblasts.obs[leiden_CAF].astype(int) + 1)
    .astype(str)
    .str.zfill(2)
)
print(f"Leiden 聚类结果已更新，并将前缀修改为 'CAF '")

sc.tl.umap(tumor_associated_fibroblasts)

umap_key = f'X_umap_CAF_n{n_neighbors}_r{resolution:.1f}'
tumor_associated_fibroblasts.obsm[umap_key] = tumor_associated_fibroblasts.obsm["X_umap"].copy()

umap_coords = tumor_associated_fibroblasts.obsm['X_umap']
tumor_associated_fibroblasts.obs[leiden_CAF] = tumor_associated_fibroblasts.obs[leiden_CAF].astype(str)

df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
df[leiden_CAF] = tumor_associated_fibroblasts.obs[leiden_CAF].values

umap_plot_filename = f'umap_leidenCAF_n{n_neighbors}_r{resolution:.1f}.png'

sc.pl.umap(
    tumor_associated_fibroblasts,
    color=[leiden_CAF],
    title=f'UMAP LeidenCAF n{n_neighbors} r{resolution:.1f}',
    save=umap_plot_filename,
    show=False,
    size=30,
    alpha=0.4,
    legend_loc='on data',
    legend_fontsize=12,
    legend_fontweight='bold'
)



rank_genes_key = leiden_CAF

if rank_genes_key not in tumor_associated_fibroblasts.obs:
    raise ValueError(f"{rank_genes_key} not found in obs. Please check the clustering step.")

sc.tl.rank_genes_groups(tumor_associated_fibroblasts, rank_genes_key, method='wilcoxon')

import re

def move_mt_rp_genes_to_end(genes):
    mt_rp_genes = [gene for gene in genes if re.match(r'^(MT|RP)', gene, re.IGNORECASE)]
    other_genes = [gene for gene in genes if not re.match(r'^(MT|RP)', gene, re.IGNORECASE)]
    return other_genes + mt_rp_genes

rank_genes_groups = tumor_associated_fibroblasts.uns['rank_genes_groups']
for group in rank_genes_groups['names'].dtype.names:
    current_genes = rank_genes_groups['names'][group]
    rank_genes_groups['names'][group] = move_mt_rp_genes_to_end(current_genes)

if 'marker_genes' in tumor_associated_fibroblasts.uns:
    marker_genes = tumor_associated_fibroblasts.uns['marker_genes']
    tumor_associated_fibroblasts.uns['marker_genes'] = move_mt_rp_genes_to_end(marker_genes)

print("Modified rank_genes_groups:")
for group in rank_genes_groups['names'].dtype.names:
    print(f"Top 5 genes for group '{group}': {rank_genes_groups['names'][group][:5]}")

if 'marker_genes' in tumor_associated_fibroblasts.uns:
    print(f"Top 5 marker genes: {tumor_associated_fibroblasts.uns['marker_genes'][:5]}")

tracksplot_save_filename = 'tracksplot_CAF.png'
sc.pl.rank_genes_groups_tracksplot(
    tumor_associated_fibroblasts,
    n_genes=3,
    save=tracksplot_save_filename,
    show=False
)

colors_key = f'{leiden_CAF}_colors'
if colors_key in tumor_associated_fibroblasts.uns:
    leiden_colors = tumor_associated_fibroblasts.uns[colors_key]
    print(f"已提取 '{colors_key}' 配色方案。")
else:
    raise KeyError(f"'{colors_key}' 未在 tumor_associated_fibroblasts.uns 中找到。")

unique_clusters = tumor_associated_fibroblasts.obs[leiden_CAF].unique()
unique_clusters_sorted = sorted(unique_clusters)

if len(leiden_colors) < len(unique_clusters_sorted):
    raise ValueError(
        f"配色方案的颜色数量少于 '{leiden_CAF}' 聚类标签的数量。请确保配色方案包含足够的颜色。"
    )

color_mapping = {cluster: color for cluster, color in zip(unique_clusters_sorted, leiden_colors)}

proportions = pd.crosstab(
    tumor_associated_fibroblasts.obs['orig.ident'],
    tumor_associated_fibroblasts.obs[leiden_CAF],
    normalize='index'
)

desired_order = ["Pr", "LN_Me", "Dis_Me"]
existing_order = [ident for ident in desired_order if ident in proportions.index]
proportions_sorted = proportions.loc[existing_order]
orig_idents_sorted = proportions_sorted.index.tolist()

plt.ioff()

fig, ax = plt.subplots(figsize=(12, 8))
bottom = pd.Series([0] * len(proportions_sorted), index=proportions_sorted.index)

for cluster in unique_clusters_sorted:
    ax.barh(
        proportions_sorted.index,
        proportions_sorted[cluster],
        left=bottom,
        color=color_mapping[cluster],
        edgecolor='black',
        label=str(cluster)
    )
    bottom += proportions_sorted[cluster]

ax.set_xlabel('Proportion', fontsize=14)
ax.set_ylabel('orig.ident', fontsize=14)
ax.set_title(f'Composition of {leiden_CAF} in Different orig.ident', fontsize=16)
ax.legend(
    title=leiden_CAF,
    bbox_to_anchor=(1.05, 1),
    loc='upper left'
)
ax.invert_yaxis()
plt.tight_layout()

output_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/F1D_CAF.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"图像已保存至：{output_path}")

adata_old = sc.read('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC/MNFHC_mEPC_CAF.h5ad')

import pandas as pd

adata_cells = adata.obs_names
adata_old_cells = adata_old.obs_names
leiden_mEPC_CAF = adata.obs['leiden_mEPC_CAF']

leiden_df = pd.DataFrame({
    'cell': adata_cells,
    'leiden_mEPC_CAF': leiden_mEPC_CAF
})

adata_old.obs = adata_old.obs.merge(leiden_df, left_index=True, right_on='cell', how='left')
print(adata_old.obs.head())

adata = adata_old
print(adata)

sc.pl.umap(
    adata,
    color='leiden_mEPC_CAF',
    title='UMAP of EPC Cells Colored by leiden_mEPC_CAF',
    save='F1F_1.png',
    show=False,
    size=30,
    legend_loc='on data',
    legend_fontsize=10,
    legend_fontweight='bold',
    add_outline=True,
    outline_width=(0.05, 0.05)
)

sc.pl.umap(
    adata,
    color='leiden_mEPC_CAF',
    title='UMAP of EPC Cells Colored by leiden_mEPC_CAF',
    save='F1F_2.png',
    show=False,
    size=30,
    legend_loc=None,
    outline_width=(0.05, 0.05)
)

for group in adata.obs['group'].unique():
    sc.pl.umap(
        adata[adata.obs['group'] == group],
        color='leiden_mEPC_CAF',
        title=f'UMAP for group {group}',
        frameon=False,
        legend_loc=None,
        size=30,
        outline_width=(0.05, 0.05),
        save=f'SF2B_{group}.png',
        show=False
    )

print(adata.obs)
adata.write('MNFHC_mEPC_CAF_GSVA.h5ad')

matplotlib.use('Agg')

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC"
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

os.chdir(work_dir)
print(f"Current working directory: {os.getcwd()}")

adata = sc.read('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/SF1E/MNFHC_allgene_Leiden_Marker_EPC_CNV.h5ad')

category_ePC_status = adata.obs.groupby(['orig.ident', 'EPC_status']).size().unstack(fill_value=0)
category_ePC_status_percentage = category_ePC_status.div(category_ePC_status.sum(axis=1), axis=0) * 100
print(category_ePC_status_percentage)

if 'EPC_status' not in adata.obs:
    raise ValueError("EPC_status column not found in the AnnData object.")
malignant_cells = adata[adata.obs['EPC_status'] == 'Malignant', :]
print(malignant_cells.obs)

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene_celltype.h5ad"
adata = sc.read(file_path)
print(adata.obs.head())

if 'cell_type' not in adata.obs:
    raise ValueError("在 AnnData 对象中未找到 'cell_type' 列。")

total_cells = adata.shape[0]
tumor_associated_fibroblasts = adata[adata.obs['cell_type'].str.contains('Cancer-associated Fibroblast', case=False, na=False), :]
fibroblast_count = tumor_associated_fibroblasts.shape[0]
non_fibroblast_count = total_cells - fibroblast_count

print(f"总细胞数: {total_cells}")
print(f"提取的含有 'fib' 的细胞数量（tumor_associated_fibroblasts）: {fibroblast_count}")
print(f"未提取的细胞数量 (非 'fib' 细胞): {non_fibroblast_count}")
print(tumor_associated_fibroblasts.obs)

combined = malignant_cells.concatenate(
    tumor_associated_fibroblasts,
    batch_key='dataset',
    batch_categories=['MalignantEPC', 'CAF']
)
print(combined)
print(combined.obs.head())

adata = combined
original_indices = adata.obs_names.tolist()

def remove_last_dash(index):
    if '-' in index:
        return '-'.join(index.split('-')[:-1])
    else:
        return index

new_indices = [remove_last_dash(idx) for idx in original_indices]
adata.obs_names = pd.Index(new_indices)

print(adata)
print(adata.obs.head())




















