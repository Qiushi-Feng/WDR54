#### ST-ssGSVA
library(Seurat)

input_dir  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata"
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (i in 1:4) {
  rds_file <- file.path(input_dir, paste0("Patient", i, "_processed.rds"))
  seurat_obj <- readRDS(rds_file)
  scale_data <- GetAssayData(object = seurat_obj, assay = "SCT", slot = "data")
  scale_mat <- as.matrix(scale_data)
  csv_file <- file.path(output_dir, paste0("Patient", i, ".csv"))
  write.csv(scale_mat, file = csv_file, quote = FALSE)
  message("已保存: ", csv_file)
}

input_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA"

for (i in 1:4) {
  fname <- sprintf("Patient%d.csv", i)
  file_path <- file.path(input_dir, fname)
  df <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  df_t <- t(df)
  prefix <- sprintf("^Patient%d_", i)
  rownames(df_t) <- sub(prefix, "", rownames(df_t))
  write.csv(df_t, file = file_path, row.names = TRUE, quote = FALSE)
  message(fname, " 已转置、去除前缀并覆写。")
}
import scanpy as sc
import pandas as pd
import os
import numpy as np
from scipy import sparse
import glob

gmt_dir   = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/MSigDB_HALLMARK"
gmt_files = glob.glob(os.path.join(gmt_dir, "*.gmt"))

def load_genes_from_gmt(path):
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                d[parts[0]] = parts[2:]
    return d

for p in range(1, 5):
    print(f"\n##### 处理 Patient{p} #####")
    h5ad   = f"/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata_python/Patient{p}.h5ad"
    csv    = f"/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA/Patient{p}.csv"
    adata  = sc.read_h5ad(h5ad)
    expr   = pd.read_csv(csv, index_col=0)
    expr   = expr.reindex(index=adata.obs_names)
    if expr.isna().values.any():
        raise ValueError(f"CSV 缺少 {expr.isna().any(1).sum()} 个细胞 – Patient{p}")
    expr   = expr.reindex(columns=adata.var_names)
    if expr.isna().values.any():
        raise ValueError(f"CSV 缺少 {expr.isna().any().sum()} 个基因 – Patient{p}")
    adata.X = sparse.csr_matrix(expr.values.astype(np.float32))
    counts  = adata.X.toarray()
    score_df = pd.DataFrame(index=adata.obs_names)
    for gmt in gmt_files:
        gene_sets = load_genes_from_gmt(gmt)
        if len(gene_sets) != 1:
            continue
        gs_name, genes = next(iter(gene_sets.items()))
        genes_in = [g for g in genes if g in adata.var_names]
        if not genes_in:
            continue
        idx = [adata.var_names.get_loc(g) for g in genes_in]
        non_zero = np.sum(counts[:, idx] > 0, axis=1)
        mask = non_zero <= 1
        adata_tmp = adata.copy()
        sc.tl.score_genes(
            adata_tmp,
            gene_list = genes_in,
            score_name = f"{gs_name}_score",
            ctrl_size  = 50,
            n_bins     = 25,
            random_state = 0,
            copy       = False,
            use_raw    = False
        )
        adata_tmp.obs.loc[mask, f"{gs_name}_score"] = np.nan
        score_df[f"{gs_name}_score"] = adata_tmp.obs[f"{gs_name}_score"]
    out_csv = f"/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA/Hallmark_score_Patient{p}.csv"
    score_df.to_csv(out_csv)
    df_clean = pd.read_csv(out_csv)
    df_clean.columns = df_clean.columns.str.replace("HALLMARK_", "").str.replace("_score", "")
    df_clean.to_csv(out_csv, index=False)
    print(f"✓ Patient{p} 处理完成，结果保存至 {out_csv}")

base_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA"

for p in range(1, 5):
    input_path  = f"{base_dir}/Hallmark_score_Patient{p}.csv"
    output_path = f"{base_dir}/Hallmark_score_Patient{p}_normalized.csv"
    df = pd.read_csv(input_path)
    for col in df.columns[1:51]:
        low, high = df[col].quantile(0.005), df[col].quantile(0.995)
        df[col] = ((df[col] - low) / (high - low)).clip(0, 1)
    df.to_csv(output_path, index=False)
    print(f"Patient{p} 归一化完成，保存为 {output_path}")

import matplotlib.pyplot as plt

base_csv   = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/SSGSVA"
base_h5ad  = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata_python"
fig_dir    = f"{base_csv}/figure"
os.makedirs(fig_dir, exist_ok=True)

routes = {
    "EPITHELIAL_MESENCHYMAL_TRANSITION": "viridis",
    "TGF_BETA_SIGNALING"               : "magma",
    "P53_PATHWAY"                      : "cividis",
    "DNA_REPAIR"                       : "plasma"
}

for i in range(1, 5):
    print(f"\n—— 处理 Patient{i} ——")
    csv_path  = f"{base_csv}/Hallmark_score_Patient{i}_normalized.csv"
    h5ad_path = f"{base_h5ad}/Patient{i}.h5ad"
    df   = pd.read_csv(csv_path, index_col=0)
    adata = sc.read_h5ad(h5ad_path)
    df = df.reindex(adata.obs_names)
    for col in df.columns[1:51]:
        adata.obs[col] = df[col].values
    print(adata.obs.iloc[:5, :5])
    for pathway, cmap in routes.items():
        fig, ax = plt.subplots(figsize=(7.7, 8))
        sc.pl.embedding(
            adata,
            basis="spatial",
            color=pathway,
            color_map=cmap,
            size=160,
            ax=ax,
            show=False
        )
        ax.invert_yaxis()
        out_png = f"{fig_dir}/Patient{i}_{pathway.split('_')[0]}.png"
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  已保存 {out_png}")
