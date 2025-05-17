## ssgsva
import scanpy as sc
import pandas as pd
import os
import harmonypy as hm
import numpy as np
import matplotlib.pyplot as plt
import shutil
import glob
import glob
import os
import pandas as pd
import numpy as np
import scanpy as sc

gene_sets_dir = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/MSigDB_HALLMARK'
gmt_files = glob.glob(os.path.join(gene_sets_dir, '*.gmt'))
if len(gmt_files) == 0:
    raise FileNotFoundError(f"在目录 {gene_sets_dir} 中未找到任何 GMT 文件。")
print(f"找到 {len(gmt_files)} 个 GMT 文件。")

score_df = pd.DataFrame(index=malignantEPC.obs_names)
counts_data = malignantEPC.X
if not isinstance(counts_data, np.ndarray):
    counts_data = counts_data.toarray()

def load_genes_from_gmt(gmt_path):
    gene_lists = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            gene_set_name = parts[0]
            genes = parts[2:]
            gene_lists[gene_set_name] = genes
    return gene_lists

for gmt_file in gmt_files:
    print(f"\n处理基因集文件: {gmt_file}")
    gene_sets = load_genes_from_gmt(gmt_file)
    if len(gene_sets) != 1:
        raise ValueError(f"GMT 文件 {gmt_file} 中包含多个基因集，请确认每个文件只包含一个基因集。")
    gene_set_name, gene_list = next(iter(gene_sets.items()))
    print(f"使用基因集: {gene_set_name} 包含 {len(gene_list)} 个基因。")
    filtered_gene_list = [gene for gene in gene_list if gene in malignantEPC.var_names]
    print(f"过滤后使用 {len(filtered_gene_list)} 个基因进行评分。")
    if len(filtered_gene_list) == 0:
        print(f"基因集 {gene_set_name} 中没有任何基因在 malignantEPC 中找到。跳过该基因集。")
        continue
    try:
        gene_indices = [malignantEPC.var_names.get_loc(gene) for gene in filtered_gene_list]
    except KeyError as e:
        print(f"基因 {e} 未在 malignantEPC.var_names 中找到，跳过该基因集。")
        continue
    non_zero_counts_per_cell = np.sum(counts_data[:, gene_indices] > 0, axis=1)
    static_threshold = 1
    print(f"静态阈值: {static_threshold} 个非零基因表达。")
    mask_less_than_or_equal_threshold = non_zero_counts_per_cell <= static_threshold
    malignantEPC_scored = malignantEPC.copy()
    sc.tl.score_genes(
        malignantEPC_scored,
        gene_list=filtered_gene_list,
        score_name=f'{gene_set_name}_score',
        ctrl_as_ref=True,
        ctrl_size=50,
        gene_pool=None,
        n_bins=25,
        random_state=0,
        copy=False,
        use_raw=False
    )
    malignantEPC_scored.obs.loc[mask_less_than_or_equal_threshold, f'{gene_set_name}_score'] = np.nan
    score_df[f'{gene_set_name}_score'] = malignantEPC_scored.obs[f'{gene_set_name}_score']
    print(f"基因集 {gene_set_name} 的评分已添加。")

output_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv'
score_df.to_csv(output_file)
print(f"\n所有基因集的评分结果已保存为 {output_file}.")

file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv'
data = pd.read_csv(file_path)
columns_to_modify = data.columns[1:51]
new_columns = [col.replace('HALLMARK_', '').replace('_score', '') for col in columns_to_modify]
data.columns.values[1:51] = new_columns
print("修改后的前几行数据：")
print(data.head())
data.to_csv(file_path, index=False)
print(f"修改后的文件已保存到: {file_path}")

gene_sets_dir = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/MSigDB_HALLMARK'
gmt_files = glob.glob(os.path.join(gene_sets_dir, '*.gmt'))
if len(gmt_files) == 0:
    raise FileNotFoundError(f"在目录 {gene_sets_dir} 中未找到任何 GMT 文件。")
print(f"找到 {len(gmt_files)} 个 GMT 文件。")

score_df = pd.DataFrame(index=caf.obs_names)
counts_data = caf.X
if not isinstance(counts_data, np.ndarray):
    counts_data = counts_data.toarray()

def load_genes_from_gmt(gmt_path):
    gene_lists = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            gene_set_name = parts[0]
            genes = parts[2:]
            gene_lists[gene_set_name] = genes
    return gene_lists

for gmt_file in gmt_files:
    print(f"\n处理基因集文件: {gmt_file}")
    gene_sets = load_genes_from_gmt(gmt_file)
    if len(gene_sets) != 1:
        raise ValueError(f"GMT 文件 {gmt_file} 中包含多个基因集，请确认每个文件只包含一个基因集。")
    gene_set_name, gene_list = next(iter(gene_sets.items()))
    print(f"使用基因集: {gene_set_name} 包含 {len(gene_list)} 个基因。")
    filtered_gene_list = [gene for gene in gene_list if gene in caf.var_names]
    print(f"过滤后使用 {len(filtered_gene_list)} 个基因进行评分。")
    if len(filtered_gene_list) == 0:
        print(f"基因集 {gene_set_name} 中没有任何基因在 caf 中找到。跳过该基因集。")
        continue
    try:
        gene_indices = [caf.var_names.get_loc(gene) for gene in filtered_gene_list]
    except KeyError as e:
        print(f"基因 {e} 未在 caf.var_names 中找到，跳过该基因集。")
        continue
    non_zero_counts_per_cell = np.sum(counts_data[:, gene_indices] > 0, axis=1)
    static_threshold = 1
    print(f"静态阈值: {static_threshold} 个非零基因表达。")
    mask_less_than_or_equal_threshold = non_zero_counts_per_cell <= static_threshold
    caf_scored = caf.copy()
    sc.tl.score_genes(
        caf_scored,
        gene_list=filtered_gene_list,
        score_name=f'{gene_set_name}_score',
        ctrl_as_ref=True,
        ctrl_size=50,
        gene_pool=None,
        n_bins=25,
        random_state=0,
        copy=False,
        use_raw=False
    )
    caf_scored.obs.loc[mask_less_than_or_equal_threshold, f'{gene_set_name}_score'] = np.nan
    score_df[f'{gene_set_name}_score'] = caf_scored.obs[f'{gene_set_name}_score']
    print(f"基因集 {gene_set_name} 的评分已添加。")

output_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv'
score_df.to_csv(output_file)
print(f"\n所有基因集的评分结果已保存为 {output_file}.")
file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv'
data = pd.read_csv(file_path)
columns_to_modify = data.columns[1:51]
new_columns = [col.replace('HALLMARK_', '').replace('_score', '') for col in columns_to_modify]
data.columns.values[1:51] = new_columns
print("修改后的前几行数据：")
print(data.head())
data.to_csv(file_path, index=False)
print(f"修改后的文件已保存到: {file_path}")


mepc_df = pd.read_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv')
caf_df = pd.read_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv')

leiden_mEPC_CAF = adata.obs['leiden_mEPC_CAF']

mepc_df['cell_name'] = mepc_df.iloc[:, 0]
caf_df['cell_name'] = caf_df.iloc[:, 0]

mepc_df = mepc_df.set_index('cell_name').join(leiden_mEPC_CAF, how='inner')
caf_df = caf_df.set_index('cell_name').join(leiden_mEPC_CAF, how='inner')

mepc_df.reset_index(inplace=True)
caf_df.reset_index(inplace=True)

mepc_df.to_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv', index=False)
caf_df.to_csv('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv', index=False)

print("mEPC 文件更新后的前五行:")
print(mepc_df.head())

print("\nCAF 文件更新后的前五行:")
print(caf_df.head())

mepc_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv'
caf_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv'

mepc_df = pd.read_csv(mepc_file)
caf_df = pd.read_csv(caf_file)

mepc_df = mepc_df.drop(mepc_df.columns[1], axis=1)
caf_df = caf_df.drop(caf_df.columns[1], axis=1)

mepc_df.to_csv(mepc_file, index=False)
caf_df.to_csv(caf_file, index=False)

print("第二列已删除并保存更新后的文件。")

input_file_mEPC = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv'
output_file_mEPC = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC_size.csv'

df_mEPC = pd.read_csv(input_file_mEPC)

cell_name_col_mEPC = 'Unnamed: 0'
cell_type_col_mEPC = 'leiden_mEPC_CAF'

score_columns_mEPC = [col for col in df_mEPC.columns if col not in [cell_name_col_mEPC, cell_type_col_mEPC]]

grouped_mEPC = df_mEPC.groupby(cell_type_col_mEPC)

result_dict_mEPC = {}
for cell_type, group in grouped_mEPC:
    valid_counts = group[score_columns_mEPC].notna().sum()
    percentages = (valid_counts / len(group)) * 100
    result_dict_mEPC[cell_type] = percentages

result_df_mEPC = pd.DataFrame(result_dict_mEPC).transpose()
result_df_mEPC = result_df_mEPC.sort_index().sort_index(axis=1)

print("\nmEPC 文件计算结果前五行:")
print(result_df_mEPC.head())

result_df_mEPC.to_csv(output_file_mEPC)
print(f"mEPC 文件的计算完成，结果已保存为 {output_file_mEPC}.")

input_file_caf = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv'
output_file_caf = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF_size.csv'

df_caf = pd.read_csv(input_file_caf)

cell_name_col_caf = 'Unnamed: 0'
cell_type_col_caf = 'leiden_mEPC_CAF'

score_columns_caf = [col for col in df_caf.columns if col not in [cell_name_col_caf, cell_type_col_caf]]

grouped_caf = df_caf.groupby(cell_type_col_caf)

result_dict_caf = {}
for cell_type, group in grouped_caf:
    valid_counts = group[score_columns_caf].notna().sum()
    percentages = (valid_counts / len(group)) * 100
    result_dict_caf[cell_type] = percentages

result_df_caf = pd.DataFrame(result_dict_caf).transpose()
result_df_caf = result_df_caf.sort_index().sort_index(axis=1)

print("\nCAF 文件计算结果前五行:")
print(result_df_caf.head())

result_df_caf.to_csv(output_file_caf)
print(f"CAF 文件的计算完成，结果已保存为 {output_file_caf}.")

input_file_caf = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF_size.csv'
input_file_mepc = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC_size.csv'

df_caf = pd.read_csv(input_file_caf)
if 'index' in df_caf.columns:
    df_caf = df_caf.drop(columns=['index'])
print("CAF 文件前五行：")
print(df_caf.head())
df_caf.to_csv(input_file_caf, index=False)

df_mepc = pd.read_csv(input_file_mepc)
if 'index' in df_mepc.columns:
    df_mepc = df_mepc.drop(columns=['index'])
print("\nmEPC 文件前五行：")
print(df_mepc.head())
df_mepc.to_csv(input_file_mepc, index=False)
print("\n文件处理完成，已覆盖源文件。")

df_caf = pd.read_csv(input_file_caf)
df_mepc = pd.read_csv(input_file_mepc)
df_caf.columns.values[0] = 'leiden_mEPC_CAF'
df_mepc.columns.values[0] = 'leiden_mEPC_CAF'
print("CAF 文件前五行：")
print(df_caf.head())
print("\nmEPC 文件前五行：")
print(df_mepc.head())
df_caf.to_csv(input_file_caf, index=False)
df_mepc.to_csv(input_file_mepc, index=False)
print("\n列名修改完成，文件已保存。")

file1 = input_file_caf
file2 = input_file_mepc

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
if not df1.columns.equals(df2.columns):
    raise ValueError("两个文件的列名不一致，无法合并。")
merged_df = pd.concat([df1, df2], ignore_index=True)
first_column = merged_df.columns[0]
merged_df_sorted = merged_df.sort_values(by=first_column, ascending=True).reset_index(drop=True)
print("Merged and Sorted DataFrame:")
print(merged_df_sorted)
output_file = "size.csv"
merged_df_sorted.to_csv(output_file, index=False)
print(f"\n合并后的数据已保存到当前目录的 '{output_file}' 文件中。")

if os.path.exists(file1):
    os.remove(file1)
    print(f"文件 {file1} 已删除。")
else:
    print(f"文件 {file1} 不存在。")
if os.path.exists(file2):
    os.remove(file2)
    print(f"文件 {file2} 已删除。")
else:
    print(f"文件 {file2} 不存在。")

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv"
output_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC_color.csv"

df = pd.read_csv(file_path)
cell_names = df.iloc[:, 0]
cell_types = df.iloc[:, 51]
gene_scores = df.iloc[:, 1:51]
sorted_cell_names = sorted(cell_names.unique())
sorted_gene_names = sorted(gene_scores.columns)

result = pd.DataFrame(index=sorted(cell_types.unique()), columns=sorted_gene_names)
for cell_type in result.index:
    type_data = gene_scores[cell_types == cell_type]

    def remove_outliers_and_average(series):
        lower_percentile = series.quantile(0.01)
        upper_percentile = series.quantile(0.99)
        filtered_data = series[(series >= lower_percentile) & (series <= upper_percentile)]
        return filtered_data.mean()

    for gene in sorted_gene_names:
        result.loc[cell_type, gene] = remove_outliers_and_average(type_data[gene])

result_normalized = result.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv"
output_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF_color.csv"

df = pd.read_csv(file_path)
cell_names = df.iloc[:, 0]
cell_types = df.iloc[:, 51]
gene_scores = df.iloc[:, 1:51]
sorted_cell_names = sorted(cell_names.unique())
sorted_gene_names = sorted(gene_scores.columns)

result = pd.DataFrame(index=sorted(cell_types.unique()), columns=sorted_gene_names)

for cell_type in result.index:
    type_data = gene_scores[cell_types == cell_type]
    def remove_outliers_and_average(series):
        lower_percentile = series.quantile(0.01)
        upper_percentile = series.quantile(0.99)
        filtered_data = series[(series >= lower_percentile) & (series <= upper_percentile)]
        return filtered_data.mean()
    for gene in sorted_gene_names:
        gene_data = type_data[gene]
        avg_score = remove_outliers_and_average(gene_data)
        result.loc[cell_type, gene] = avg_score

result_normalized = result.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)
result_normalized.to_csv(output_file)

file_path_mEPC = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC_color.csv"
file_path_CAF = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF_color.csv"

df_mEPC = pd.read_csv(file_path_mEPC)
df_CAF = pd.read_csv(file_path_CAF)

df_mEPC.columns.values[0] = 'leiden_mEPC_CAF'
df_CAF.columns.values[0] = 'leiden_mEPC_CAF'

df_mEPC.to_csv(file_path_mEPC, index=False)
df_CAF.to_csv(file_path_CAF, index=False)

output_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/color.csv"

df_mEPC = pd.read_csv(file_path_mEPC)
df_CAF = pd.read_csv(file_path_CAF)

df_combined = pd.concat([df_mEPC, df_CAF], axis=0, ignore_index=True)
df_combined_sorted = df_combined.sort_values(by=df_combined.columns[0])
df_combined_sorted.to_csv(output_file, index=False)

os.remove(file_path_mEPC)
os.remove(file_path_CAF)

base_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/"
color_file = os.path.join(base_path, "color.csv")
size_file = os.path.join(base_path, "size.csv")
output_plot = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/F1E_1.png"

color_data = pd.read_csv(color_file, index_col=0)
size_data = pd.read_csv(size_file, index_col=0)

color_long = color_data.reset_index().melt(id_vars=color_data.index.name, var_name='Gene', value_name='Color')
size_long = size_data.reset_index().melt(id_vars=size_data.index.name, var_name='Gene', value_name='Size')

color_long.rename(columns={color_data.index.name: 'Sample'}, inplace=True)
size_long.rename(columns={size_data.index.name: 'Sample'}, inplace=True)

merged_data = pd.merge(color_long, size_long, on=['Sample', 'Gene'])

plt.figure(figsize=(20, 10))
sns.scatterplot(
    data=merged_data,
    x='Gene',
    y='Sample',
    hue='Color',
    size='Size',
    palette='viridis',
    sizes=(20, 200),
    alpha=1.0,
    edgecolor='black',
    linewidth=0.5
)
plt.xticks(rotation=90)
plt.xlabel('Gene')
plt.ylabel('Sample')
plt.title('Dotplot of Gene Expression')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Legend', borderaxespad=0.)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"Dotplot 已成功保存到 {output_plot}")
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/figures/F1E_2.png"

fig, ax = plt.subplots(figsize=(6, 1))
fig.subplots_adjust(bottom=0.5)

cmap = mpl.cm.get_cmap('viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=1)

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, orientation='horizontal', ax=ax)
cbar.set_label('Viridis Color Scale')

ax.set_axis_off()

plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Viridis 色条已成功保存到 {output_path}")

import pandas as pd

caf_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_CAF.csv'
mepc_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/Hallmark_score_mEPC.csv'

caf_df = pd.read_csv(caf_file)
mepc_df = pd.read_csv(mepc_file)

combined_df = pd.concat([caf_df, mepc_df], axis=0, ignore_index=True)

combined_df = combined_df.drop(columns=['leiden_mEPC_CAF'])

print("Columns after merging and removing 'leiden_mEPC_CAF':")
print(combined_df.columns)

output_file = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/File8_GSVA_for_mEPC_and_CAF.csv'
combined_df.to_csv(output_file, index=False)

print(f"File saved to {output_file}")

## GO analysis
import scanpy as sc
import pandas as pd

input_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC/MNFHC_CAF.h5ad"
output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2/CAF_DEG.csv"
adata = sc.read(input_path)
if 'rank_genes_groups' not in adata.uns:
    raise KeyError("'rank_genes_groups' 不存在于 adata.uns 中。请确保已运行相关的差异表达分析。")
rank_genes_groups = adata.uns['rank_genes_groups']
groups = rank_genes_groups['names'].dtype.names
deg_list = []
for group in groups:
    genes = rank_genes_groups['names'][group]
    logfc = rank_genes_groups['logfoldchanges'][group]
    pvals = rank_genes_groups['pvals'][group]
    pvals_adj = rank_genes_groups['pvals_adj'][group]
    scores = rank_genes_groups['scores'][group]
    temp_df = pd.DataFrame({
        'cell_type': group,
        'gene': genes,
        'logfoldchange': logfc,
        'pval': pvals,
        'pval_adj': pvals_adj,
        'score': scores
    })
    deg_list.append(temp_df)
deg_df = pd.concat(deg_list, ignore_index=True)
deg_df.to_csv(output_path, index=False)
print(f"差异表达基因数据已成功保存至：{output_path}")

input_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC/MNFHC_MalignantEPC.h5ad"
output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2/mEPC_DEG.csv"
adata = sc.read(input_path)
if 'rank_genes_groups' not in adata.uns:
    raise KeyError("'rank_genes_groups' 不存在于 adata.uns 中。请确保已运行相关的差异表达分析。")
rank_genes_groups = adata.uns['rank_genes_groups']
groups = rank_genes_groups['names'].dtype.names
deg_list = []
for group in groups:
    genes = rank_genes_groups['names'][group]
    logfc = rank_genes_groups['logfoldchanges'][group]
    pvals = rank_genes_groups['pvals'][group]
    pvals_adj = rank_genes_groups['pvals_adj'][group]
    scores = rank_genes_groups['scores'][group]
    temp_df = pd.DataFrame({
        'cell_type': group,
        'gene': genes,
        'logfoldchange': logfc,
        'pval': pvals,
        'pval_adj': pvals_adj,
        'score': scores
    })
    deg_list.append(temp_df)
deg_df = pd.concat(deg_list, ignore_index=True)
deg_df.to_csv(output_path, index=False)
print(f"差异表达基因数据已成功保存至：{output_path}")

setwd("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2")
library("org.Hs.eg.db")
rt = read.csv("CAF_DEG.csv", header = TRUE, check.names = FALSE)
genes = as.vector(rt$gene)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt$entrezID <- entrezIDs
write.csv(rt, file = "CAF_DEG_annotation.csv", row.names = FALSE, quote = FALSE)

rt = read.csv("mEPC_DEG.csv", header = TRUE, check.names = FALSE)
genes = as.vector(rt$gene)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt$entrezID <- entrezIDs
write.csv(rt, file = "mEPC_DEG_annotation.csv", row.names = FALSE, quote = FALSE)

file.remove("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2/CAF_DEG.csv")
file.remove("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2/mEPC_DEG.csv")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF2")
if (!dir.exists("figures")) {
    dir.create("figures")
}
rt <- read.csv("CAF_DEG_annotation.csv", header = TRUE, sep = ",", check.names = FALSE)
print(colnames(rt))
print(head(rt, 5))

selected_genes <- rt[rt$pval < 0.01 & rt$logfoldchange > 0.5 & !is.na(rt$entrezID) & rt$cell_type == 'CAF 04', ]
gene_list <- as.character(selected_genes$entrezID)

kk <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  ont = "all",
  readable = TRUE
)

figures_path <- "figures"

png(filename = file.path(figures_path, "SF2C_CAF 04.png"), width = 1000, height = 800, res = 150)
barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

write.table(kk, file = file.path("GO_CAF 04.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

selected_genes <- rt[rt$pval < 0.01 & rt$logfoldchange > 0.5 & !is.na(rt$entrezID) & rt$cell_type == 'CAF 06', ]
gene_list <- as.character(selected_genes$entrezID)

kk <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  ont = "all",
  readable = TRUE
)

write.table(kk, file = "GO_Selected_Genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

figures_path <- "figures"

png(filename = file.path(figures_path, "SF2C_CAF 06.png"), width = 1000, height = 800, res = 150)
barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

write.table(kk, file = file.path("GO_CAF 06.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!dir.exists("figures")) {
  dir.create("figures")
}

rt <- read.csv("mEPC_DEG_annotation.csv", header = TRUE, sep = ",", check.names = FALSE)
print(colnames(rt))
print(head(rt, 5))
selected_genes <- rt[rt$pval < 0.01 & rt$logfoldchange > 0.5 & !is.na(rt$entrezID) & rt$cell_type == 'mEPC 08', ]

gene_list <- as.character(selected_genes$entrezID)

kk <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  ont = "all",
  readable = TRUE
)

figures_path <- "figures"

png(filename = file.path(figures_path, "SF2C_mEPC 08.png"), width = 1000, height = 800, res = 150)
barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

write.table(kk, file = file.path("GO_mEPC 08.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
