### cellchat
import os
import scanpy as sc
import pandas as pd

work_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1"
os.chdir(work_dir)

print("当前工作路径:", os.getcwd())

input_file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF1/MNFHC_Leiden_allgene.h5ad"
adata = sc.read(input_file)
adata.obs = adata.obs.drop('is_outlier', axis=1)

new_cluster_names = [
    'T Cell 01', 'Macrophage', 'Natural Killer Cell', 'Epithelial Cell 01', 'Monocyte',
    'Plasma Cell 01', 'Epithelial Cell 02', 'Cancer-associated Fibroblast 01', 'Epithelial Cell 03',
    'B cell', 'Endothelial Cell', 'T Cell 02', 'Cancer-associated Fibroblast 02', 'Epithelial Cell 04',
    'T Cell 03', 'Cancer-associated Fibroblast 03', 'Plasmacytoid Dendritic Cell', 'Epithelial Cell 05',
    'Salivary cell', 'Mast Cell', 'Langerhans Cell', 'CD8-positive, Alpha-Beta T Cell',
    'Plasma Cell 02', 'Smooth Muscle Cell', 'Mononuclear Macrophage', 'Ciliated Cell'
]

adata.rename_categories('leiden_n10_r0.4', new_cluster_names)

print("leiden_n10_r0.4 与细胞类型对应关系:")
print(adata.obs[['leiden_n10_r0.4']].value_counts())

df = pd.DataFrame({
    "cell": adata.obs.index,
    "leiden_n10_r0.4": adata.obs["leiden_n10_r0.4"].astype(str)
})

out_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "cellclusters.csv")

df.to_csv(out_path, index=False)

print(f"已保存细胞聚类信息至: {out_path}")
adata = sc.read_h5ad("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1EF/MNFHC_mEPC_CAF_GSVA.h5ad")

df = adata.obs[['leiden_mEPC_CAF']].reset_index()
df.columns = ['cell', 'leiden_mEPC_CAF']

output_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/leiden_mEPC_CAF.csv"
df.to_csv(output_path, index=False)

print(df.head())

csv_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/leiden_mEPC_CAF.csv"
df = pd.read_csv(csv_path)

df['cell'] = (
    df['cell']
    .str.replace('LN_', '', regex=False)
    .str.replace('Dis_', '', regex=False)
)

df.to_csv(csv_path, index=False)

print(df.head())
library(Seurat)
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F0/scRNA/Processed/Processed_Rawcounts/step2_harmony/HarmonyResults_OldSeuratforAnndata.rds")

csv_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/cellclusters.csv"
cluster_df <- read.csv(csv_path, stringsAsFactors = FALSE)
rownames(cluster_df) <- cluster_df$cell

seurat_obj@meta.data$leiden_clusters <- cluster_df[rownames(seurat_obj@meta.data), "leiden_n10_r0.4"]

head(seurat_obj@meta.data, 5)

out_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/MNFHC_immunecellchat.rds"
saveRDS(seurat_obj, file = out_path)
message("已保存 Seurat 对象至: ", out_path)

mapping <- read.csv(
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/leiden_mEPC_CAF.csv",
  stringsAsFactors = FALSE
)

meta <- seurat_obj@meta.data
meta$cellchat_rough <- NA_character_

idx_map <- match(rownames(meta), mapping$cell)
has_map <- !is.na(idx_map)
meta$cellchat_rough[has_map] <- mapping$leiden_mEPC_CAF[idx_map[has_map]]

lc <- meta$leiden_clusters
not_filled <- is.na(meta$cellchat_rough)

cond_tcell <- not_filled & (grepl("T Cell 0", lc) | lc == "CD8-positive, Alpha-Beta T Cell")
meta$cellchat_rough[cond_tcell] <- "T Cell"

cond_bcell <- not_filled & lc == "B cell"
meta$cellchat_rough[cond_bcell] <- "B cell"

cond_plasma <- not_filled & grepl("Plasma Cell", lc)
meta$cellchat_rough[cond_plasma] <- "Plasma Cell"

cond_mono_macro <- not_filled & lc %in% c("Mononuclear Macrophage", "Macrophage")
meta$cellchat_rough[cond_mono_macro] <- "Mononuclear Macrophage"

other_types <- c(
  "Plasmacytoid Dendritic Cell",
  "Monocyte",
  "Natural Killer Cell",
  "Mast Cell",
  "Langerhans Cell"
)
cond_other <- not_filled & lc %in% other_types
meta$cellchat_rough[cond_other] <- lc[cond_other]

meta$cellchat_rough[is.na(meta$cellchat_rough)] <- "None"

seurat_obj@meta.data <- meta
saveRDS(
  seurat_obj,
  file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/MNFHC_immunecellchat.rds"
)

n_before <- nrow(seurat_obj@meta.data)
message("删除前细胞数: ", n_before)

keep_cells <- rownames(seurat_obj@meta.data)[
  seurat_obj@meta.data$cellchat_rough != "None" &
  !grepl("CAF", seurat_obj@meta.data$cellchat_rough)
]

seurat_filtered <- subset(seurat_obj, cells = keep_cells)

n_after <- nrow(seurat_filtered@meta.data)
message("删除后细胞数: ", n_after)

saveRDS(
  seurat_filtered,
  file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/MNFHC_cellchatrough.rds"
)

seurat_obj <- seurat_filtered
rm(seurat_filtered)

figures_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/rough"
if (!dir.exists(figures_path)) {
    dir.create(figures_path, recursive = TRUE)
}
setwd(figures_path)

library(patchwork)
library(viridis)
options(Seurat.object.assay.version = "v3")
library(CellChat)

data.input <- seurat_obj@assays$RNA@data
meta <- seurat_obj@meta.data
cell.use <- rownames(meta)
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$cellchat_rough)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_rough")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellchat_rough")

levels(cellchat@idents)

groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

future::plan("multicore", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/rough/cellchat_rough.rds"
saveRDS(cellchat, file = save_path)
cat("CellChat object has been saved to", save_path, "\n")

color_mapping <- c(
  "mEPC 01"                     = "#E377C2",
  "mEPC 02"                     = "#B5BD61",
  "mEPC 03"                     = "#17BEDF",
  "mEPC 04"                     = "#AEC7E8",
  "mEPC 05"                     = "#FFBB78",
  "mEPC 06"                     = "#98DF8A",
  "mEPC 07"                     = "#FF9896",
  "mEPC 08"                     = "#C5B0D5",
  "mEPC 09"                     = "#C49C94",
  "mEPC 10"                     = "#F7B6D2",
  "T Cell"                      = "#BC8F8F",
  "B cell"                      = "#E5AF44",
  "Plasma Cell"                 = "#AFE544",
  "Plasmacytoid Dendritic Cell"= "#44E544",
  "Monocyte"                    = "#44E5AF",
  "Mononuclear Macrophage"      = "#44AFE5",
  "Natural Killer Cell"         = "#D3D3D3",
  "Langerhans Cell"             = "#AF44E5",
  "Mast Cell"                   = "#E544AF"
)

groups_in_network <- levels(cellchat@idents)
color.use <- color_mapping[groups_in_network]

if (any(is.na(color.use))) {
  warning("Some group labels have no matching color: ",
          paste(groups_in_network[is.na(color.use)], collapse = ", "))
}

png("F4_mEPC_vs_Others.png", width = 5000, height = 5000, res = 300)

net <- cellchat@net$weight
group_names <- rownames(net)

has_mEPC <- grepl("mEPC", group_names)

mask <- matrix(FALSE,
               nrow = length(group_names),
               ncol = length(group_names),
               dimnames = list(group_names, group_names))

mask[has_mEPC, !has_mEPC] <- TRUE
mask[!has_mEPC, has_mEPC] <- TRUE

filtered_net <- net
filtered_net[!mask] <- 0

netVisual_circle(
  filtered_net,
  vertex.weight    = groupSize,
  weight.scale     = TRUE,
  label.edge       = FALSE,
  edge.width.max   = 15,
  arrow.width      = 2,
  arrow.size       = 0.5,
  color.use        = color.use,
  vertex.label.cex = 0.5
)

dev.off()

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/rough/seperation"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

net <- cellchat@net$weight
group_names <- rownames(net)

mEPC_groups <- grep("mEPC", group_names, value = TRUE)

for (grp in mEPC_groups) {
  filename <- file.path(out_dir, paste0(grp, ".png"))
  png(
    filename = filename,
    width    = 5000,
    height   = 5000,
    units    = "px",
    res      = 300
  )
  
  is_grp  <- group_names == grp
  non_grp <- !grepl("mEPC", group_names)
  mask <- matrix(
    FALSE,
    nrow     = length(group_names),
    ncol     = length(group_names),
    dimnames = list(group_names, group_names)
  )
  mask[is_grp,  non_grp] <- TRUE
  mask[non_grp, is_grp ] <- TRUE
  
  filtered_net <- net
  filtered_net[!mask] <- 0
  
  netVisual_circle(
    filtered_net,
    vertex.weight    = groupSize,
    weight.scale     = TRUE,
    label.edge       = FALSE,
    edge.width.max   = 15,
    arrow.width      = 2,
    arrow.size       = 0.5,
    color.use        = color.use,
    vertex.label.cex = 0.5
  )
  
  dev.off()
}

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/rough"
png(
  filename = file.path(out_dir, "net_heatmap_count.png"),
  width    = 2000,
  height   = 2000,
  units    = "px",
  res      = 300
)
gg_count <- netVisual_heatmap(
  cellchat,
  measure    = "count",
  color.use  = color_mapping,
  color.heatmap = c("#0D0887FF","#F0F921FF")
)
print(gg_count)
dev.off()

png(
  filename = file.path(out_dir, "net_heatmap_weight.png"),
  width    = 2000,
  height   = 2000,
  units    = "px",
  res      = 300
)
gg_weight <- netVisual_heatmap(
  cellchat,
  measure    = "weight",
  color.use  = color_mapping,
  color.heatmap = c("#0D0887FF","#F0F921FF")
)
print(gg_weight)
dev.off()

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/rough/ligand_receptor"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

groups        <- levels(cellchat@idents)
mEPC_groups   <- grep("mEPC", groups, value = TRUE)
non_mEPC_grps <- setdiff(groups, mEPC_groups)

for (src in non_mEPC_grps) {
  out_file <- file.path(out_dir, paste0(gsub(" ", "_", src), "_to_mEPC.png"))
  png(filename = out_file, width = 1800, height = 1800, units = "px", res = 300)
  
  p <- netVisual_bubble(
    object      = cellchat,
    sources.use = src,
    targets.use = mEPC_groups,
    angle.x     = 45
  )
  print(p)
  dev.off()
  message("Saved bubble plot for ", src, " -> mEPC groups to:\n  ", out_file)
}
groups <- levels(cellchat@idents)
src <- "mEPC 08"
targets <- setdiff(groups, grep("mEPC", groups, value = TRUE))
out_file <- file.path(out_dir, paste0(gsub(" ", "_", src), "_to_non_mEPC.png"))
png(filename = out_file, width = 1900, height = 1900, units = "px", res = 300)
p <- netVisual_bubble(
  object      = cellchat,
  sources.use = src,
  targets.use = targets,
  angle.x     = 45
)
print(p)
dev.off()
message("Saved bubble plot for ", src, " → non-mEPC groups to:\n  ", out_file)

out_dir <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/rough/ligand_receptor"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
groups <- levels(cellchat@idents)
sources <- setdiff(groups, grep("mEPC", groups, value = TRUE))
target <- "mEPC 08"
out_file <- file.path(out_dir, paste0("non_mEPC_to_", gsub(" ", "_", target), ".png"))
png(filename = out_file, width = 1900, height = 1900, units = "px", res = 300)
p <- netVisual_bubble(
  object      = cellchat,
  sources.use = sources,
  targets.use = target,
  angle.x     = 45
)
print(p)
dev.off()
message("Saved bubble plot for non-mEPC → ", target, " to:\n  ", out_file)

cells_to_keep <- WhichCells(seurat_obj, expression = cellchat_rough %in% c("mEPC 08", "Mononuclear Macrophage"))
seurat_detailed <- subset(seurat_obj, cells = cells_to_keep)
saveRDS(
  seurat_detailed,
  file = "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/MNFHC_cellchatdetailed.rds"
)

library(Seurat)
library(harmony)

seurat_ob <- readRDS("/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/MNFHC_cellchatdetailed.rds")
cells_to_keep <- rownames(seurat_ob@meta.data)[!grepl("mEPC", seurat_ob@meta.data$cellchat_rough)]
seurat_obj <- subset(seurat_ob, cells = cells_to_keep)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE)
seurat_obj <- RunPCA(
  seurat_obj,
  assay = "RNA",
  features = VariableFeatures(seurat_obj),
  npcs = 50,
  verbose = FALSE
)
cat("[INFO] PCA 完成。\n")
seurat_obj <- RunHarmony(
  object         = seurat_obj,
  group.by.vars  = "batch",
  assay.use      = "RNA",
  reduction       = "pca",
  reduction.save  = "harmony",
  dims.use        = 1:50,
  verbose         = FALSE
)
cat("[INFO] Harmony 去批次完成\n")
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims      = 1:30,
  k.param   = 100,
  verbose   = FALSE
)
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.8,
  algorithm  = 4,
  verbose    = FALSE
)
cat("[INFO] 邻域构建与聚类完成。\n")
print(table(seurat_obj@meta.data$seurat_clusters))
detailed_clusters <- data.frame(
  cell            = rownames(seurat_obj@meta.data),
  seurat_clusters = seurat_obj@meta.data$seurat_clusters,
  row.names       = NULL,
  stringsAsFactors= FALSE
)
output_file <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/detailed_clusters.csv"
write.csv(detailed_clusters, file = output_file, row.names = FALSE)
cat("[INFO] 已将细胞索引及其 seurat_clusters 保存至：", output_file, "\n")

library(Seurat)
library(dplyr)
all_markers <- FindAllMarkers(
  object           = seurat_obj,
  only.pos         = TRUE,
  min.pct          = 0.1,
  logfc.threshold  = 0.25
)
top100_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  slice_head(n = 100) %>%
  ungroup()
output_file <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/detailed/top100_DEGs_per_cluster.csv"
write.csv(top100_markers, file = output_file, row.names = FALSE)
cat("[INFO] 已将每个簇前100个差异表达基因保存至：", output_file, "\n")

library(dplyr)
file_path <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/detailed_clusters.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)
cluster_map <- c(
  "1"  = "Mononuclear Macrophage",
  "2"  = "Classical Macrophage",
  "3"  = "Myeloid Progenitor Cell",
  "4"  = "Classical Macrophage",
  "5"  = "Fibroblast",
  "6"  = "Mononuclear Macrophage",
  "7"  = "Epithelial Cell",
  "8"  = "Mononuclear Macrophage",
  "9"  = "Classical Monocyte",
  "10" = "Mature NK T cell",
  "11" = "Myeloid Progenitor Cell",
  "12" = "Classical Monocyte",
  "13" = "Fibroblast"
)
df <- df %>%
  mutate(
    detailed_clusters = cluster_map[as.character(seurat_clusters)]
  )
write.csv(df, file = file_path, row.names = FALSE)
head(df, 5)

library(Seurat)
library(dplyr)
rds_path <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/MNFHC_cellchatdetailed.rds"
seurat_ob <- readRDS(rds_path)
csv_path <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/detailed_clusters.csv"
df_clusters <- read.csv(csv_path, stringsAsFactors = FALSE)
cluster_map <- setNames(df_clusters$detailed_clusters, df_clusters$cell)
meta <- seurat_ob@meta.data
meta$cell <- rownames(meta)
meta <- meta %>%
  mutate(
    cellchat_detailed = case_when(
      cellchat_rough == "mEPC 08"     ~ "mEPC 08",
      cell %in% names(cluster_map)    ~ cluster_map[cell],
      TRUE                             ~ "none"
    )
  )
meta$cell <- NULL
seurat_ob@meta.data <- meta
saveRDS(seurat_ob, file = rds_path)
cat("[INFO] 已成功写入 cellchat_detailed 列并覆盖保存至：", rds_path, "\n")

seurat_full <- readRDS("/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/MNFHC_cellchatdetailed.rds")
to_remove <- c("Mature NK T cell", "Epithelial Cell", "Fibroblast")
keep_cells <- rownames(seurat_full@meta.data)[
  ! seurat_full@meta.data$cellchat_detailed %in% to_remove
]
seurat_obj <- subset(seurat_full, cells = keep_cells)
print(table(seurat_obj@meta.data$cellchat_detailed))

figures_path <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/detailed"
if (!dir.exists(figures_path)) dir.create(figures_path, recursive = TRUE)
setwd(figures_path)
library(patchwork)
library(Seurat)
library(viridis)
options(Seurat.object.assay.version = "v3")
library(CellChat)
data.input <- seurat_obj@assays$RNA@.data
meta <- seurat_obj@meta.data
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$cellchat_detailed)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_detailed")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellchat_detailed")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save_path <- "/media/desk16/tly6105/Essential_and_single_cell/F4/cellchat/rough/cellchat_rough.rds"
saveRDS(cellchat, file = save_path)
saveRDS(cellchat, file = save_path)
cat("CellChat object has been saved to", save_path, "\n")

color_mapping <- c(
  "mEPC 08"                 = "#C5B0D5",
  "Classical Macrophage"    = "#AA3A49",
  "Myeloid Progenitor Cell" = "#E19D49",
  "Mononuclear Macrophage"  = "#4E8872",
  "Classical Monocyte"      = "#9F8DB8"
)
groups_in_network <- levels(cellchat@idents)
color.use <- color_mapping[groups_in_network]
if (any(is.na(color.use))) {
  warning("有些群体标签在 color_mapping 中没有对应的颜色：",
          paste(groups_in_network[is.na(color.use)], collapse = ", "),
          "\n请检查映射向量是否完整或拼写一致。")
}

png("F4_mEPC 08_vs_Others.png", width = 5000, height = 5000, res = 300)
net <- cellchat@net$weight
group_names <- rownames(net)
has_mEPC <- grepl("mEPC", group_names)
mask <- matrix(FALSE,
               nrow = length(group_names),
               ncol = length(group_names),
               dimnames = list(group_names, group_names))
mask[has_mEPC, !has_mEPC] <- TRUE
mask[!has_mEPC, has_mEPC] <- TRUE
filtered_net <- net
filtered_net[!mask] <- 0
netVisual_circle(
  filtered_net,
  vertex.weight    = groupSize,
  weight.scale     = TRUE,
  label.edge       = FALSE,
  edge.width.max   = 15,
  arrow.width      = 2,
  arrow.size       = 0.5,
  color.use        = color.use,
  vertex.label.cex = 0.5
)
dev.off()

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/cellchat/detailed"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
groups <- levels(cellchat@idents)
mEPC_groups <- grep("mEPC", groups, value = TRUE)
non_mEPC_grps <- setdiff(groups, mEPC_groups)
for (src in non_mEPC_grps) {
  out_file <- file.path(out_dir, paste0(gsub(" ", "_", src), "_to_mEPC 08.png"))
  png(filename = out_file, width = 1800, height = 1800, units = "px", res = 300)
  p <- netVisual_bubble(
    object      = cellchat,
    sources.use = src,
    targets.use = mEPC_groups,
    angle.x     = 45
  )
  print(p)
  dev.off()
  message("Saved bubble plot for ", src, " -> mEPC groups to:\n  ", out_file)
}

groups <- levels(cellchat@idents)
src <- "mEPC 08"
targets <- setdiff(groups, grep("mEPC", groups, value = TRUE))
out_file <- file.path(out_dir, paste0(gsub(" ", "_", src), "_to_non_mEPC.png"))
png(filename = out_file, width = 1900, height = 1900, units = "px", res = 300)
p <- netVisual_bubble(
  object      = cellchat,
  sources.use = src,
  targets.use = targets,
  angle.x     = 45
)
print(p)
dev.off()
