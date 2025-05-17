### Anndata construction
## R part
library(Seurat)
library(Matrix)
input_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts"
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
for(i in 1:4){
  patient_id <- paste0("Patient", i)
  input_file <- file.path(input_dir, paste0(patient_id, "_processed.rds"))
  cat("正在处理：", patient_id, "\n")
  obj <- readRDS(input_file)
  DefaultAssay(obj) <- "SCT"
  counts_mat <- GetAssayData(obj, slot = "counts", assay = "SCT")
  mtx_file <- file.path(output_dir, paste0(patient_id, "_matrix.mtx"))
  writeMM(counts_mat, file = mtx_file)
  barcodes_file <- file.path(output_dir, paste0(patient_id, "_barcodes.tsv"))
  write.table(colnames(counts_mat), file = barcodes_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  features_file <- file.path(output_dir, paste0(patient_id, "_features.tsv"))
  write.table(rownames(counts_mat), file = features_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  clusters_file <- file.path(output_dir, paste0(patient_id, "_clusters.tsv"))
  meta <- obj@meta.data
  cluster_df <- data.frame(barcode = rownames(meta), seurat_clusters = meta$seurat_clusters, stringsAsFactors = FALSE)
  write.table(cluster_df, file = clusters_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat(patient_id, "处理完成。\n")
}
cat("所有患者的 SCT raw counts 文件均已保存至：", output_dir, "\n")
## python part
import os
import json
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.io import mmread
from PIL import Image
base_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts"
file_suffixes = ["barcodes.tsv", "clusters.tsv"]
for patient_num in range(1, 5):
    patient_prefix = f"Patient{patient_num}_"
    for suffix in file_suffixes:
        file_path = os.path.join(base_dir, patient_prefix + suffix)
        if os.path.exists(file_path):
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            new_content = content.replace(patient_prefix, "")
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            print(f"已处理文件：{file_path}")
        else:
            print(f"文件不存在：{file_path}")
print("所有文件中的相应前缀已删除并成功覆盖！")
mtx_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts/Patient1_matrix.mtx"
barcodes_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts/Patient1_barcodes.tsv"
features_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts/Patient1_features.tsv"
image_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/image/Patient1_image"
positions_path = os.path.join(image_dir, "tissue_positions_list.csv")
scale_factors_path = os.path.join(image_dir, "scalefactors_json.json")
save_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient1.h5ad"
print("Reading matrix file...")
mtx = mmread(mtx_path).T.tocsr()
print("Reading barcodes file...")
barcodes = pd.read_csv(barcodes_path, header=None)
barcodes.columns = ["barcode"]
print("Reading features file...")
features = pd.read_csv(features_path, header=None)
features.columns = ["gene_id"]
print("Constructing AnnData...")
adata = anndata.AnnData(X=mtx)
adata.obs["barcode"] = barcodes["barcode"].values
adata.obs_names = adata.obs["barcode"]
adata.var["gene_id"] = features["gene_id"].values
adata.var_names = adata.var["gene_id"]
print("Reading spatial positions...")
pos = pd.read_csv(positions_path, header=None)
pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
pos.set_index("barcode", inplace=True)
pos = pos.loc[adata.obs_names, :]
adata.obsm["spatial"] = pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
adata.obs["in_tissue"] = pos["in_tissue"].values
adata.obs["array_row"] = pos["array_row"].values
adata.obs["array_col"] = pos["array_col"].values
print("Reading scale factors...")
with open(scale_factors_path, "r") as f:
    scale_factors = json.load(f)
library_id = "Patient1"
adata.uns["spatial"] = {library_id: {"scalefactors": scale_factors, "images": {}}}
image_files = {"detected_tissue_image": "detected_tissue_image.jpg", "tissue_hires_image": "tissue_hires_image.png", "aligned_fiducials": "aligned_fiducials.jpg"}
for key, filename in image_files.items():
    img_path = os.path.join(image_dir, filename)
    if os.path.exists(img_path):
        with Image.open(img_path) as img:
            adata.uns["spatial"][library_id]["images"][key] = np.array(img)
    else:
        print(f"Warning: {filename} not found, skipping.")
print("Writing AnnData to h5ad...")
adata.write_h5ad(save_path)
print("Finished! Your AnnData is saved to:", save_path)
import anndata
import pandas as pd
h5ad_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient1.h5ad"
clusters_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/SCT_rawcounts/Patient1_clusters.tsv"
adata = anndata.read_h5ad(h5ad_path)
clusters_df = pd.read_csv(clusters_path, sep="\t")
clusters_df.set_index("barcode", inplace=True)
clusters_df["seurat_clusters"] = clusters_df["seurat_clusters"].astype(str)
common_cells = adata.obs.index.intersection(clusters_df.index)
adata.obs.loc[common_cells, "seurat_clusters"] = clusters_df.loc[common_cells, "seurat_clusters"]
print(f"更新了 {len(common_cells)} 个细胞的 seurat_clusters 信息。")
adata.write_h5ad(h5ad_path)
print("seurat_clusters 信息已保存到：", h5ad_path)
import os
import scanpy as sc
data_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python"
patient_mappings = {
    "Patient1.h5ad": {"1": "#2171A8", "2": "#EC7C1E", "3": "#37955F", "4": "#642D84", "5": "#D74692", "6": "#EADE25", "7": "#8C564B", "8": "#6AC3CA", "9": "#C5B0D5", "10": "#F4B579", "11": "#98CC85", "12": "#D62728", "13": "#B2DF8A", "14": "#BABABA"},
    "Patient2.h5ad": {"1": "#2171A8", "2": "#EC7C1E", "3": "#37955F", "4": "#642D84", "5": "#C5B0D5", "6": "#D74692", "7": "#EADE25", "8": "#D62728", "9": "#6AC3CA", "10": "#F4B579", "11": "#8C564B", "12": "#98CC85", "13": "#B2DF8A", "14": "#BABABA", "15": "#2F2D54", "16": "#E0A4DD"},
    "Patient3.h5ad": {"1": "#2171A8", "2": "#EC7C1E", "3": "#D62728", "4": "#C5B0D5", "5": "#8C564B", "6": "#37955F", "7": "#642D84", "8": "#D74692", "9": "#EADE25", "10": "#6AC3CA", "11": "#F4B579", "12": "#98CC85"},
    "Patient4.h5ad": {"1": "#2171A8", "2": "#C5B0D5", "3": "#EC7C1E", "4": "#37955F", "5": "#642D84", "6": "#D62728", "7": "#8C564B", "8": "#D74692", "9": "#EADE25", "10": "#6AC3CA", "11": "#F4B579", "12": "#98CC85", "13": "#B2DF8A"}
}

for fn, cmap in patient_mappings.items():
    fp = os.path.join(data_dir, fn)
    adata = sc.read_h5ad(fp)
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype(str)
    order = (adata.obs["seurat_clusters"].cat.categories.tolist() if pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]) else sorted(adata.obs["seurat_clusters"].unique(), key=lambda x: int(x)))
    adata.uns["seurat_clusters_colors"] = [cmap.get(cl, "#CCCCCC") for cl in order]
    adata.write(fp)

### stlearn
import scanpy as sc
import stlearn as st
import os
import matplotlib.pyplot as plt
file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient1.h5ad'
data = sc.read_h5ad(file_path)
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn")
data = st.convert_scanpy(data)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.raw = data
st.pp.scale(data)
st.em.run_pca(data, n_comps=50, random_state=0)
library_id = "Patient1"
if "hires" not in data.uns["spatial"][library_id]["images"]:
    data.uns["spatial"][library_id]["images"]["hires"] = data.uns["spatial"][library_id]["images"]["tissue_hires_image"]
st.pp.tiling(data, out_path="Patient1_test", crop_size=100)
st.pp.extract_feature(data)
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=30, use_rep='X_pca_morphology', random_state=0)
out_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn/Patient1"
data.uns["iroot"] = st.spatial.trajectory.set_root(
    data,
    use_label="seurat_clusters",
    cluster="2",
    use_raw=True
)
st.spatial.trajectory.pseudotime(
    data,
    eps=50,
    use_rep="X_pca",
    use_label="seurat_clusters"
)
st.spatial.trajectory.pseudotimespace_global(
    data,
    use_label="seurat_clusters",
    list_clusters=["9", "7", "12"]
)
st.pl.cluster_plot(
    data,
    use_label="seurat_clusters",
    show_trajectories=True,
    list_clusters=["9", "7", "12"],
    show_subcluster=False,
    size=25,
    figsize=(8, 8),
    color_bar_size=15,
    cell_alpha=0.85,
    trajectory_width=1.5,
    trajectory_arrowsize=9,
    reverse=False,
    show_cluster_labels=True,
    trajectory_edge_color="black",
    show_plot=False,
    crop=False,
    show_color_bar=False,
    text_box_size=12
)
plot1_path = os.path.join(out_dir, "plot1.png")
plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
plt.close()
st.spatial.trajectory.pseudotimespace_global(
    data,
    use_label="seurat_clusters",
    list_clusters=["9", "7"]
)
st.pl.cluster_plot(
    data,
    use_label="seurat_clusters",
    show_trajectories=True,
    list_clusters=["9", "7"],
    show_subcluster=False,
    size=25,
    figsize=(8, 8),
    color_bar_size=15,
    cell_alpha=0.85,
    reverse=False,
    crop=False,
    show_cluster_labels=True,
    trajectory_edge_color="black",
    show_plot=False,
    show_color_bar=False,
    text_box_size=12
)
plot2_path = os.path.join(out_dir, "plot2.png")
plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
plt.close()
st.spatial.trajectory.pseudotimespace_global(
    data,
    use_label="seurat_clusters",
    list_clusters=["12", "7"]
)
st.pl.cluster_plot(
    data,
    use_label="seurat_clusters",
    show_trajectories=True,
    list_clusters=["12", "7"],
    show_subcluster=False,
    size=25,
    figsize=(8, 8),
    color_bar_size=15,
    cell_alpha=0.85,
    reverse=False,
    crop=False,
    show_cluster_labels=True,
    trajectory_edge_color="black",
    show_plot=False,
    show_color_bar=False,
    text_box_size=12
)
plot3_path = os.path.join(out_dir, "plot3.png")
plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
plt.close()
st.spatial.trajectory.pseudotimespace_global(
    data,
    use_label="seurat_clusters",
    list_clusters=["12", "9"]
)
st.pl.cluster_plot(
    data,
    use_label="seurat_clusters",
    show_trajectories=True,
    list_clusters=["12", "9"],
    show_subcluster=False,
    size=25,
    figsize=(8, 8),
    color_bar_size=15,
    cell_alpha=0.85,
    reverse=False,
    crop=False,
    show_cluster_labels=True,
    trajectory_edge_color="black",
    show_plot=False,
    show_color_bar=False,
    text_box_size=12
)
plot4_path = os.path.join(out_dir, "plot4.png")
plt.savefig(plot4_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"已保存 iroot = {iroot_val} 对应的图像于文件夹：{out_dir}")
file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient2.h5ad'
data = sc.read_h5ad(file_path)
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn")
data = st.convert_scanpy(data)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.raw = data
st.pp.scale(data)
st.em.run_pca(data, n_comps=50, random_state=0)
library_id = "Patient2"
if "hires" not in data.uns["spatial"][library_id]["images"]:
    data.uns["spatial"][library_id]["images"]["hires"] = data.uns["spatial"][library_id]["images"]["tissue_hires_image"]
st.pp.tiling(data, out_path="Patient2_test", crop_size=100)
st.pp.extract_feature(data)
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=30, use_rep='X_pca_morphology', random_state=0)
output_root = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn/Patient2"
for root_val in range(1, 17):
    folder_name = f"iroot_{root_val}"
    out_dir = os.path.join(output_root, folder_name)
    os.makedirs(out_dir, exist_ok=True)
    data.uns["iroot"] = st.spatial.trajectory.set_root(
        data,
        use_label="seurat_clusters",
        cluster=str(root_val),
        use_raw=True
    )
    st.spatial.trajectory.pseudotime(
        data,
        eps=50,
        use_rep="X_pca",
        use_label="seurat_clusters"
    )
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["8", "5", "11"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["8", "5", "11"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        trajectory_width=1.5,
        trajectory_arrowsize=9,
        reverse=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        crop=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot1_path = os.path.join(out_dir, f"iroot_{root_val}_plot1.png")
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["8", "5"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["8", "5"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot2_path = os.path.join(out_dir, f"iroot_{root_val}_plot2.png")
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["11", "5"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["11", "5"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot3_path = os.path.join(out_dir, f"iroot_{root_val}_plot3.png")
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["11", "8"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["11", "8"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot4_path = os.path.join(out_dir, f"iroot_{root_val}_plot4.png")
    plt.savefig(plot4_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"已保存 iroot = {root_val} 对应的图像于文件夹：{out_dir}")
print("所有起点图像均已保存！")
file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient3.h5ad'
data = sc.read_h5ad(file_path)
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn")
data = st.convert_scanpy(data)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.raw = data
st.pp.scale(data)
st.em.run_pca(data, n_comps=50, random_state=0)
library_id = "Patient3"
if "hires" not in data.uns["spatial"][library_id]["images"]:
    data.uns["spatial"][library_id]["images"]["hires"] = data.uns["spatial"][library_id]["images"]["tissue_hires_image"]
st.pp.tiling(data, out_path="Patient3_test", crop_size=100)
st.pp.extract_feature(data)
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=30, use_rep='X_pca_morphology', random_state=0)
output_root = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn/Patient3"
for root_val in range(1, 13):
    folder_name = f"iroot_{root_val}"
    out_dir = os.path.join(output_root, folder_name)
    os.makedirs(out_dir, exist_ok=True)
    data.uns["iroot"] = st.spatial.trajectory.set_root(
        data,
        use_label="seurat_clusters",
        cluster=str(root_val),
        use_raw=True
    )
    st.spatial.trajectory.pseudotime(
        data,
        eps=50,
        use_rep="X_pca",
        use_label="seurat_clusters"
    )
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["4", "5", "3"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["4", "5", "3"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        trajectory_width=1.5,
        trajectory_arrowsize=9,
        reverse=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        crop=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot1_path = os.path.join(out_dir, f"iroot_{root_val}_plot1.png")
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["4", "5"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["4", "5"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot2_path = os.path.join(out_dir, f"iroot_{root_val}_plot2.png")
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["3", "5"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["3", "5"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot3_path = os.path.join(out_dir, f"iroot_{root_val}_plot3.png")
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["3", "4"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["3", "4"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot4_path = os.path.join(out_dir, f"iroot_{root_val}_plot4.png")
    plt.savefig(plot4_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"已保存 iroot = {root_val} 对应的图像于文件夹：{out_dir}")
print("所有起点图像均已保存！")
file_path = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient4.h5ad'
data = sc.read_h5ad(file_path)
os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn")
data = st.convert_scanpy(data)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.raw = data
st.pp.scale(data)
st.em.run_pca(data, n_comps=50, random_state=0)
library_id = "Patient4"
if "hires" not in data.uns["spatial"][library_id]["images"]:
    data.uns["spatial"][library_id]["images"]["hires"] = data.uns["spatial"][library_id]["images"]["tissue_hires_image"]
st.pp.tiling(data, out_path="Patient4_test", crop_size=100)
st.pp.extract_feature(data)
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=30, use_rep='X_pca_morphology', random_state=0)
output_root = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stlearn/Patient4"
for root_val in range(1, 14):
    folder_name = f"iroot_{root_val}"
    out_dir = os.path.join(output_root, folder_name)
    os.makedirs(out_dir, exist_ok=True)
    data.uns["iroot"] = st.spatial.trajectory.set_root(
        data,
        use_label="seurat_clusters",
        cluster=str(root_val),
        use_raw=True
    )
    st.spatial.trajectory.pseudotime(
        data,
        eps=50,
        use_rep="X_pca",
        use_label="seurat_clusters"
    )
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["7", "6", "2"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["7", "6", "2"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        trajectory_width=2.0,
        trajectory_arrowsize=12,
        reverse=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        crop=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot1_path = os.path.join(out_dir, f"iroot_{root_val}_plot1.png")
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["7", "6"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["7", "6"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot2_path = os.path.join(out_dir, f"iroot_{root_val}_plot2.png")
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["2", "6"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["2", "6"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot3_path = os.path.join(out_dir, f"iroot_{root_val}_plot3.png")
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    plt.close()
    st.spatial.trajectory.pseudotimespace_global(
        data,
        use_label="seurat_clusters",
        list_clusters=["2", "7"]
    )
    st.pl.cluster_plot(
        data,
        use_label="seurat_clusters",
        show_trajectories=True,
        list_clusters=["2", "7"],
        show_subcluster=False,
        size=25,
        figsize=(8, 8),
        color_bar_size=15,
        cell_alpha=0.85,
        reverse=False,
        crop=False,
        show_cluster_labels=True,
        trajectory_edge_color="black",
        show_plot=False,
        show_color_bar=False,
        text_box_size=12
    )
    plot4_path = os.path.join(out_dir, f"iroot_{root_val}_plot4.png")
    plt.savefig(plot4_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"已保存 iroot = {root_val} 对应的图像于文件夹：{out_dir}")
print("所有起点图像均已保存！")

### SPAtrack
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import scanpy as sc
import os
import json
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.io import mmread
from PIL import Image
import spaTrack as spt

sc.settings.verbosity = 0
plt.rcParams["figure.dpi"] = 300

h5ad_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient1.h5ad"
adata = anndata.read_h5ad(h5ad_path)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.obs["cell_type"] = adata.obs["seurat_clusters"]
adata.obs["cluster"] = adata.obs["seurat_clusters"]

start_cells = spt.set_start_cells(
    adata,
    select_way='cell_type',
    cell_type="9",
    basis='spatial',
    split=False,
    n_neigh=5
)
print("使用 set_start_cells 选中的起始细胞索引：", start_cells)

adata.obsm["X_spatial"] = adata.obsm["spatial"].astype(np.float64).copy()
adata.obsm["X_spatial"][:, 0] *= 0.8

adata.obsp['trans'] = spt.get_ot_matrix(adata, data_type='spatial', alpha1=0, alpha2=1)
adata.obs['ptime'] = spt.get_ptime(adata, start_cells)
adata.uns['E_grid'], adata.uns['V_grid'] = spt.get_velocity(adata, basis='spatial', n_neigh_pos=100)

fig1, ax1 = plt.subplots(figsize=(7.7, 8))
sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='ptime',
    show=False,
    ax=ax1,
    color_map='Reds',
    title='Ptime',
    size=130
)
ax1.invert_yaxis()
output_path1 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient1_ptime.png"
plt.savefig(output_path1, dpi=300, bbox_inches='tight')
plt.close(fig1)
print("ptime 图已保存到:", output_path1)
fig2, ax2 = plt.subplots(figsize=(7.2, 8))


adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')


colors_list = list(adata.uns["seurat_clusters_colors"])


cats = adata.obs['seurat_clusters'].cat.categories
palette_dict = dict(zip(cats, colors_list))


sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='seurat_clusters',
    palette=palette_dict,
    show=False,
    ax=ax2,
    legend_loc=None,
    frameon=False,
    title='Trajectory',
    alpha=0.8,
    size=130
)

# 叠加流线图
ax2.streamplot(
    adata.uns['E_grid'][0],
    adata.uns['E_grid'][1],
    adata.uns['V_grid'][0],
    adata.uns['V_grid'][1],
    color='black',
    linewidth=2,
    density=1.2,
    arrowsize=1.2
)
ax2.invert_yaxis()
output_path2 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient1_trajectory.png"
plt.savefig(output_path2, dpi=300, bbox_inches='tight')
plt.close(fig2)
print("Trajectory 图已保存到:", output_path2)

sc.settings.verbosity = 0
plt.rcParams["figure.dpi"] = 300

h5ad_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient2.h5ad"
adata = anndata.read_h5ad(h5ad_path)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.obs["cell_type"] = adata.obs["seurat_clusters"]
adata.obs["cluster"] = adata.obs["seurat_clusters"]

start_cells = spt.set_start_cells(
    adata,
    select_way='cell_type',
    cell_type="5",
    basis='spatial',
    split=False,
    n_neigh=5
)
print("使用 set_start_cells 选中的起始细胞索引：", start_cells)

adata.obsm["X_spatial"] = adata.obsm["spatial"].astype(np.float64).copy()
adata.obsm["X_spatial"][:, 0] *= 0.8

adata.obsp['trans'] = spt.get_ot_matrix(adata, data_type='spatial', alpha1=0, alpha2=1)
adata.obs['ptime'] = spt.get_ptime(adata, start_cells)
adata.uns['E_grid'], adata.uns['V_grid'] = spt.get_velocity(adata, basis='spatial', n_neigh_pos=100)

fig1, ax1 = plt.subplots(figsize=(7.7, 8))
sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='ptime',
    show=False,
    ax=ax1,
    color_map='Reds',
    title='Ptime',
    size=130
)
ax1.invert_yaxis()
output_path1 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient2_ptime.png"
plt.savefig(output_path1, dpi=300, bbox_inches='tight')
plt.close(fig1)
print("ptime 图已保存到:", output_path1)
fig2, ax2 = plt.subplots(figsize=(7.2, 8))


adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')


colors_list = list(adata.uns["seurat_clusters_colors"])


cats = adata.obs['seurat_clusters'].cat.categories
palette_dict = dict(zip(cats, colors_list))


sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='seurat_clusters',
    palette=palette_dict,
    show=False,
    ax=ax2,
    legend_loc=None,
    frameon=False,
    title='Trajectory',
    alpha=0.8,
    size=130
)

# 叠加流线图
ax2.streamplot(
    adata.uns['E_grid'][0],
    adata.uns['E_grid'][1],
    adata.uns['V_grid'][0],
    adata.uns['V_grid'][1],
    color='black',
    linewidth=2,
    density=1.2,
    arrowsize=1.2
)
ax2.invert_yaxis()
output_path2 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient2_trajectory.png"
plt.savefig(output_path2, dpi=300, bbox_inches='tight')
plt.close(fig2)
print("Trajectory 图已保存到:", output_path2)

sc.settings.verbosity = 0
plt.rcParams["figure.dpi"] = 300

h5ad_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient3.h5ad"
adata = anndata.read_h5ad(h5ad_path)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.obs["cell_type"] = adata.obs["seurat_clusters"]
adata.obs["cluster"] = adata.obs["seurat_clusters"]

start_cells = spt.set_start_cells(
    adata,
    select_way='cell_type',
    cell_type="4",
    basis='spatial',
    split=False,
    n_neigh=5
)
print("使用 set_start_cells 选中的起始细胞索引：", start_cells)

adata.obsm["X_spatial"] = adata.obsm["spatial"].astype(np.float64).copy()
adata.obsm["X_spatial"][:, 0] *= 0.8

adata.obsp['trans'] = spt.get_ot_matrix(adata, data_type='spatial', alpha1=0, alpha2=1)
adata.obs['ptime'] = spt.get_ptime(adata, start_cells)
adata.uns['E_grid'], adata.uns['V_grid'] = spt.get_velocity(adata, basis='spatial', n_neigh_pos=100)

fig1, ax1 = plt.subplots(figsize=(7.7, 8))
sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='ptime',
    show=False,
    ax=ax1,
    color_map='Reds',
    title='Ptime',
    size=130
)
ax1.invert_yaxis()
output_path1 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient3_ptime.png"
plt.savefig(output_path1, dpi=300, bbox_inches='tight')
plt.close(fig1)
print("ptime 图已保存到:", output_path1)
fig2, ax2 = plt.subplots(figsize=(7.2, 8))


adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')


colors_list = list(adata.uns["seurat_clusters_colors"])


cats = adata.obs['seurat_clusters'].cat.categories
palette_dict = dict(zip(cats, colors_list))


sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='seurat_clusters',
    palette=palette_dict,
    show=False,
    ax=ax2,
    legend_loc=None,
    frameon=False,
    title='Trajectory',
    alpha=0.8,
    size=130
)

# 叠加流线图
ax2.streamplot(
    adata.uns['E_grid'][0],
    adata.uns['E_grid'][1],
    adata.uns['V_grid'][0],
    adata.uns['V_grid'][1],
    color='black',
    linewidth=2,
    density=1.2,
    arrowsize=1.2
)
ax2.invert_yaxis()
output_path2 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient3_trajectory.png"
plt.savefig(output_path2, dpi=300, bbox_inches='tight')
plt.close(fig2)
print("Trajectory 图已保存到:", output_path2)

sc.settings.verbosity = 0
plt.rcParams["figure.dpi"] = 300

h5ad_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient4.h5ad"
adata = anndata.read_h5ad(h5ad_path)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.obs["cell_type"] = adata.obs["seurat_clusters"]
adata.obs["cluster"] = adata.obs["seurat_clusters"]

start_cells = spt.set_start_cells(
    adata,
    select_way='cell_type',
    cell_type="2",
    basis='spatial',
    split=False,
    n_neigh=5
)
print("使用 set_start_cells 选中的起始细胞索引：", start_cells)

adata.obsm["X_spatial"] = adata.obsm["spatial"].astype(np.float64).copy()
adata.obsm["X_spatial"][:, 0] *= 0.8

adata.obsp['trans'] = spt.get_ot_matrix(adata, data_type='spatial', alpha1=0, alpha2=1)
adata.obs['ptime'] = spt.get_ptime(adata, start_cells)
adata.uns['E_grid'], adata.uns['V_grid'] = spt.get_velocity(adata, basis='spatial', n_neigh_pos=100)

fig1, ax1 = plt.subplots(figsize=(7.7, 8))
sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='ptime',
    show=False,
    ax=ax1,
    color_map='Reds',
    title='Ptime',
    size=130
)
ax1.invert_yaxis()
output_path1 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient4_ptime.png"
plt.savefig(output_path1, dpi=300, bbox_inches='tight')
plt.close(fig1)
print("ptime 图已保存到:", output_path1)
fig2, ax2 = plt.subplots(figsize=(7.2, 8))


adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')


colors_list = list(adata.uns["seurat_clusters_colors"])


cats = adata.obs['seurat_clusters'].cat.categories
palette_dict = dict(zip(cats, colors_list))


sc.pl.embedding(
    adata,
    basis='X_spatial',
    color='seurat_clusters',
    palette=palette_dict,
    show=False,
    ax=ax2,
    legend_loc=None,
    frameon=False,
    title='Trajectory',
    alpha=0.8,
    size=130
)

# 叠加流线图
ax2.streamplot(
    adata.uns['E_grid'][0],
    adata.uns['E_grid'][1],
    adata.uns['V_grid'][0],
    adata.uns['V_grid'][1],
    color='black',
    linewidth=2,
    density=1.2,
    arrowsize=1.2
)
ax2.invert_yaxis()
output_path2 = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/spatrack/Patient4_trajectory.png"
plt.savefig(output_path2, dpi=300, bbox_inches='tight')
plt.close(fig2)
print("Trajectory 图已保存到:", output_path2)