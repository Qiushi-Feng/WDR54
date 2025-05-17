### ST-monocle
library(Seurat)

base_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata"
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle2"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

mapping_list <- list(
  "1" = c("12" = "CAF 04", "7" = "CAF 06", "9" = "mEPC 08"),
  "2" = c("8"  = "CAF 04", "11" = "CAF 06", "5" = "mEPC 08"),
  "3" = c("3"  = "CAF 04", "5" = "CAF 06", "4"  = "mEPC 08"),
  "4" = c("6"  = "CAF 04", "7" = "CAF 06", "2"  = "mEPC 08")
)

subset_list <- list()

for (i in 1:4) {
  file_path <- file.path(base_path, paste0("Patient", i, "_processed.rds"))
  seurat_obj <- readRDS(file_path)
  current_mapping <- mapping_list[[as.character(i)]]
  allowed_clusters <- names(current_mapping)
  subset_obj <- subset(seurat_obj, subset = seurat_clusters %in% allowed_clusters)
  subset_obj@meta.data$leiden <- current_mapping[as.character(subset_obj@meta.data$seurat_clusters)]
  subset_obj@meta.data$Patient <- paste0("Patient", i)
  subset_list[[i]] <- subset_obj
  cat("Patient", i, "处理完成，共筛选出", ncol(subset_obj), "个细胞。\n")
}

merged_obj <- subset_list[[1]]
if (length(subset_list) > 1) {
  for (j in 2:length(subset_list)) {
    merged_obj <- merge(merged_obj, y = subset_list[[j]])
  }
}

output_file <- file.path(output_dir, "merged_for_monocle.rds")
saveRDS(merged_obj, file = output_file)
cat("合并后的 Seurat 对象已保存到:", output_file, "\n")

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle2/merged_for_monocle.rds")
outputDir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/monocle2"

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

cat("Output directory is set to:", outputDir, "\n")
set.seed(123)
DefaultAssay(seurat_obj) <- 'SCT'
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)

gbm_cds <- newCellDataSet(
  as(as.matrix(seurat_obj@assays[["SCT"]]@counts), 'sparseMatrix'),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
expressed_genes <- row.names(subset(fData(gbm_cds)))
disp_table <- dispersionTable(gbm_cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1.0 * dispersion_fit)
ordering_genes <- ordering_genes_temp$gene_id
gbm_cds <- setOrderingFilter(gbm_cds, ordering_genes)
write.table(ordering_genes, file = paste0(outputDir, "/", "all_dispersion.ordering_genes.tsv"), 
            row.names = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
png(file = paste(outputDir, "all_ordering_genes.png", sep = '/'), width = 800, height = 600)
print(plot_ordering_genes(gbm_cds))
dev.off()

gbm_cds <- reduceDimension(
  gbm_cds,
  max_components = 2,
  method = "DDRTree"
)

gbm_cds <- orderCells(gbm_cds)
gbm_cds$barcode <- colnames(gbm_cds)
rep <- as.data.frame(pData(gbm_cds))
write.table(pData(gbm_cds), file = paste0(outputDir, "/", "cell_Pseudotime.txt"), row.names = T, quote = F, sep = ',')
state_levels <- levels(pData(gbm_cds)$State)
if (length(state_levels) <= 2) {  
  widths_state = 8  
  heights_state = 5  
  nrows_state = 1
} else if (1 < length(state_levels) & length(state_levels) <= 8) {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 2
} else {  
  widths_state = length(state_levels) * 1.5  
  heights_state = length(state_levels) * 1.5  
  nrows_state = 3
}

celltype_levels <- levels(gbm_cds$celltype)
if (length(celltype_levels) <= 2) {  
  widths_celltype = 8  
  heights_celltype = 5  
  nrows_celltype = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 2
} else {  
  widths_celltype = length(celltype_levels) * 1.5  
  heights_celltype = length(celltype_levels) * 1.5  
  nrows_celltype = 3
}

sample_levels <- levels(gbm_cds$orig.ident)
if (length(celltype_levels) <= 2) {  
  widths_orig = 8  
  heights_orig = 5  
  nrows_orig = 1
} else if (1 < length(celltype_levels) & length(celltype_levels) <= 8) {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 2
} else {  
  widths_orig = length(celltype_levels) * 1.5  
  heights_orig = length(celltype_levels) * 1.5  
  nrows_orig = 3
}

p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
png(file = paste0(outputDir, "/1_Pseudotime_trajectory.png"), width = 800, height = 600)
print(p1)
dev.off()

p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
png(file = paste0(outputDir, "/1_Pseudotime_trajectory_legend.png"), width = 800, height = 600)
print(p1)
dev.off()

p1 <- plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = 3) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
      facet_wrap(~leiden)
png(file = paste0(outputDir, "/2_Pseudotime_trajectory_splited_by_orig.png"), width = 1400, height = 1200)
print(p1)
dev.off()

MYCOLOR <- c("#D62728", "#8C564B", "#C5B0D5")
p3 <- plot_cell_trajectory(gbm_cds, color_by = 'leiden', cell_size = 3) + 
  scale_color_manual(values = scales::alpha(MYCOLOR, 0.6)) + 
  ggtitle('celltype') + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
png(file = paste0(outputDir, "/4_Celltype_trajectory.png"), width = 800, height = 600)
print(p3)
dev.off()

### ST-stavia
import scanpy as sc
import omicverse as ov
from omicverse.externel import VIA
import matplotlib.pyplot as plt
ov.plot_set()
import os
import pandas as pd
import numpy as np
np.int = int
import harmonypy as hm

os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stavia")
print("当前工作目录为：", os.getcwd())

adata = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient1.h5ad')
adata.X = adata.X.astype('float64')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)

hvgs = adata.var_names[adata.var['highly_variable_features']].tolist()

adata = adata[:, hvgs]
ov.pp.scale(adata)

ncomps = 30
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_pca'
clusters = 'seurat_clusters'

v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['9']
)
v0.run_VIA()

color_mapping = {
    "1":  "#2171A8",
    "2":  "#EC7C1E",
    "3":  "#37955F",
    "4":  "#642D84",
    "5":  "#D74692",
    "6":  "#EADE25",
    "7":  "#8C564B",
    "8":  "#6AC3CA",
    "9":  "#C5B0D5",
    "10": "#F4B579",
    "11": "#98CC85",
    "12": "#D62728",
    "13": "#B2DF8A",
    "14": "#BABABA"
}

if not pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]):
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype("category")

cluster_order = list(adata.obs["seurat_clusters"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["seurat_clusters_colors"] = color_list

fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="seurat_clusters",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)

output_path = os.path.join(os.getcwd(), "Patient1_viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

print(f"图形已保存至: {output_path}")




os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stavia")
print("当前工作目录为：", os.getcwd())

adata = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient2.h5ad')
adata.X = adata.X.astype('float64')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)

hvgs = adata.var_names[adata.var['highly_variable_features']].tolist()

adata = adata[:, hvgs]
ov.pp.scale(adata)

ncomps = 30
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_pca'
clusters = 'seurat_clusters'

v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['5']
)
v0.run_VIA()

color_mapping = {
    "1":  "#2171A8",
    "2":  "#EC7C1E",
    "3":  "#37955F",
    "4":  "#642D84",
    "5":  "#C5B0D5",
    "6":  "#D74692",
    "7":  "#EADE25",
    "8":  "#D62728",
    "9":  "#6AC3CA",
    "10": "#F4B579",
    "11": "#8C564B",
    "12": "#98CC85",
    "13": "#B2DF8A",
    "14": "#BABABA",
    "15": "#2F2D54",
    "16": "#E0A4DD"
}


if not pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]):
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype("category")

cluster_order = list(adata.obs["seurat_clusters"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["seurat_clusters_colors"] = color_list

fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="seurat_clusters",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)

output_path = os.path.join(os.getcwd(), "Patient2_viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

print(f"图形已保存至: {output_path}")


os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stavia")
print("当前工作目录为：", os.getcwd())

adata = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient3.h5ad')
adata.X = adata.X.astype('float64')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)

hvgs = adata.var_names[adata.var['highly_variable_features']].tolist()

adata = adata[:, hvgs]
ov.pp.scale(adata)

ncomps = 30
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_pca'
clusters = 'seurat_clusters'

v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['4']
)
v0.run_VIA()

color_mapping = {
    "1":  "#2171A8",
    "2":  "#EC7C1E",
    "3":  "#D62728",
    "4":  "#C5B0D5",
    "5":  "#8C564B",
    "6":  "#37955F",
    "7":  "#642D84",
    "8":  "#D74692",
    "9":  "#EADE25",
    "10": "#6AC3CA",
    "11": "#F4B579",
    "12": "#98CC85"
}



if not pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]):
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype("category")

cluster_order = list(adata.obs["seurat_clusters"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["seurat_clusters_colors"] = color_list

fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="seurat_clusters",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)

output_path = os.path.join(os.getcwd(), "Patient3_viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

print(f"图形已保存至: {output_path}")


os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stavia")
print("当前工作目录为：", os.getcwd())

adata = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient4.h5ad')
adata.X = adata.X.astype('float64')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)

hvgs = adata.var_names[adata.var['highly_variable_features']].tolist()

adata = adata[:, hvgs]
ov.pp.scale(adata)

ncomps = 30
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_pca'
clusters = 'seurat_clusters'

v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['4']
)
v0.run_VIA()


os.chdir("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/stavia")
print("当前工作目录为：", os.getcwd())

adata = sc.read_h5ad('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata_python/Patient4.h5ad')
adata.X = adata.X.astype('float64')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50)

adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)

hvgs = adata.var_names[adata.var['highly_variable_features']].tolist()

adata = adata[:, hvgs]
ov.pp.scale(adata)

ncomps = 30
knn = 15
v0_random_seed = 4
memory = 10
dataset = ''
use_rep = 'X_pca'
clusters = 'seurat_clusters'

v0 = VIA.core.VIA(
    data=adata.obsm[use_rep][:, 0:ncomps],
    true_label=adata.obs[clusters],
    edgepruning_clustering_resolution=0.15,
    cluster_graph_pruning=0.15,
    knn=knn,
    resolution_parameter=1.5,
    dataset=dataset,
    random_seed=v0_random_seed,
    memory=memory,
    root_user=['2']
)
v0.run_VIA()

color_mapping = {
    "1":  "#2171A8",
    "2":  "#C5B0D5",
    "3":  "#EC7C1E",
    "4":  "#37955F",
    "5":  "#642D84",
    "6":  "#D62728",
    "7":  "#8C564B",
    "8":  "#D74692",
    "9":  "#EADE25",
    "10": "#6AC3CA",
    "11": "#F4B579",
    "12": "#98CC85",
    "13": "#B2DF8A"
}




if not pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]):
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype("category")

cluster_order = list(adata.obs["seurat_clusters"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["seurat_clusters_colors"] = color_list

fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="seurat_clusters",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)

output_path = os.path.join(os.getcwd(), "Patient4_viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

print(f"图形已保存至: {output_path}")




if not pd.api.types.is_categorical_dtype(adata.obs["seurat_clusters"]):
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype("category")

cluster_order = list(adata.obs["seurat_clusters"].cat.categories)
color_list = [color_mapping.get(cl, "#000000") for cl in cluster_order]
adata.uns["seurat_clusters_colors"] = color_list

fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(
    adata,
    clusters="seurat_clusters",
    dpi=300,
    via_object=v0,
    ax_text=False,
    show_legend=False
)
fig.set_size_inches(16, 8)

output_path = os.path.join(os.getcwd(), "Patient4_viagraph_piechart.png")
fig.savefig(output_path, format='png', dpi=300)
plt.close(fig)

print(f"图形已保存至: {output_path}")
