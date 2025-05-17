## pyscenic
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import loompy as lp
import subprocess
import sys
import os
import scipy.sparse as sp
import scipy
from adjustText import adjust_text

sc.settings.verbosity = 3
sc.logging.print_header()
sc.set_figure_params(dpi=600, dpi_save=600)

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
from pyscenic.cli.utils import load_signatures

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2C/MNFHC_witch_counts.h5ad"
adata = sc.read_h5ad(file_path)
raw_X = adata.raw.X

if sp.issparse(raw_X):
    adata.X = raw_X.copy()
else:
    adata.X = raw_X.copy()

print("已使用 adata.raw.X 覆盖 adata.X，新的 adata.X 形状为：", adata.X.shape)

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=25)

f_loom_path_scenic = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/MNFHC_countsloom_pyscenic.loom"
hs_tfs_fname = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/reference/allTFs_hg38.txt"
adjacencies_fname = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/MNFHC_ptscenic_adjacencies.TSV"
reference_folder = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/reference"
filenames = [
    "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
    "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
]
ranking_dbs_fnames = [os.path.join(reference_folder, fn) for fn in filenames]
dbs_param = ' '.join(ranking_dbs_fnames)
motif_annotations_fname = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/reference/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
motifs_fname = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/MNFHC_motif.csv"
outloom = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/MNFHC_scenicresult.loom"
aucell_mtx_fname = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/MNFHC_scenic_aucell.csv"

row_attrs = {
    "Gene": np.array(adata.var_names)
}
col_attrs = {
    "CellID": np.array(adata.obs_names),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten()
}

lp.create(f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)

cmd = [
    "pyscenic", "grn",
    f_loom_path_scenic,
    hs_tfs_fname,
    "-o", adjacencies_fname,
    "--num_workers", "16"
]
proc = subprocess.run(cmd, capture_output=True, text=True)
print(proc.stdout)
if proc.returncode != 0:
    print(proc.stderr, file=sys.stderr)
    sys.exit(proc.returncode)

subprocess.run(
    f"pyscenic ctx {adjacencies_fname} {dbs_param} "
    f"--annotations_fname {motif_annotations_fname} "
    f"--expression_mtx_fname {f_loom_path_scenic} "
    f"--output {motifs_fname} "
    "--num_workers 16",
    shell=True,
    check=True
)

subprocess.run(
    f"pyscenic aucell {f_loom_path_scenic} {motifs_fname} "
    f"--output {outloom} --num_workers 16",
    shell=True,
    check=True
)

subprocess.run([
    "pyscenic", "aucell",
    f_loom_path_scenic,
    motifs_fname,
    "--output", aucell_mtx_fname,
    "--num_workers", "16"
], check=True)

auc_mtx = pd.read_csv(aucell_mtx_fname, index_col=0)

regulons = load_signatures(motifs_fname)

adata = sc.read('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2C/MNFHC_witch_counts.h5ad')
add_scenic_metadata(adata, auc_mtx, regulons)
df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['leiden_mEPC_CAF']]
df_results = (
    (df_scores.groupby(by='leiden_mEPC_CAF').mean() - df_obs[signature_column_names].mean())
    / df_obs[signature_column_names].std()
).stack().reset_index().rename(columns={'level_1': 'regulon', 0: 'Z'})
df_results['Regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))

output_dir = '/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic'
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'MNFHC_regulons_results.csv')
df_results[df_results.Z >= 1].sort_values(['leiden_mEPC_CAF', 'Z'], ascending=False).to_csv(output_file, index=False)

df_heatmap = pd.pivot_table(
    data=df_results[df_results.Z >= 1].sort_values('Z', ascending=False),
    index='leiden_mEPC_CAF',
    columns='Regulon',
    values='Z'
)
fig, ax = plt.subplots(1, 1, figsize=(24, 10))
sns.heatmap(
    df_heatmap,
    ax=ax,
    annot=True,
    fmt=".1f",
    linewidths=.7,
    cbar=False,
    square=True,
    linecolor='gray',
    vmin=1,
    vmax=3,
    cmap="Blues",
    annot_kws={"size": 8}
)
ax.set_ylabel('Z score')

output_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/figure"
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "heatmap_regulons_z_score.png")
fig.tight_layout()
fig.savefig(output_file, dpi=600, format='png')

rss = regulon_specificity_scores(auc_mtx, adata.obs['leiden_mEPC_CAF'])
output_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/figure"
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "MNFHC_SCENIC_RSS.csv")
rss.T.to_csv(output_file)

cats = sorted(list(set(adata.obs['leiden_mEPC_CAF'])))
sc.set_figure_params(dpi=600, dpi_save=600)

import os
import matplotlib.pyplot as plt
from adjustText import adjust_text

output_dir = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/pyscenic/figure"
os.makedirs(output_dir, exist_ok=True)

cats = [
    'CAF 04', 'CAF 06', 'mEPC 08']
dpi_val = 600

for c in cats:
    fig, ax = plt.subplots(figsize=(6, 5))
    x = rss.T[c]
    plot_rss(rss, c, top_n=30, max_n=None, ax=ax)
    margin = (x.max() - x.min()) * 0.05
    ax.set_ylim(x.min() - margin, x.max() + margin)
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_xlabel('')
    ax.set_ylabel('Regulon specificity score (RSS)')
    ax.set_title(c)
    adjust_text(
        ax.texts,
        autoalign='xy',
        ha='right',
        va='bottom',
        arrowprops=dict(arrowstyle='-', color='lightgrey'),
        precision=0.001
    )
    fname = f"plots_regulons_{c.replace(' ', '_')}.png"
    output_path = os.path.join(output_dir, fname)
    fig.savefig(output_path, dpi=dpi_val, format='png')
    plt.close(fig)
    print(f"Figure saved to: {output_path}")

## cellchat
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(CellChat)

figures_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/figures"
if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}
setwd(figures_path)

seurat_obj <- readRDS('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/MNFHC_for_cellchat.rds')
data.input <- seurat_obj@assays$RNA@data
meta <- seurat_obj@meta.data
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "leiden_mEPC_CAF")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "leiden_mEPC_CAF")

groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
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

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/cellchat_all_signal.rds"
saveRDS(cellchat, file = save_path)
cat("CellChat object has been saved to", save_path, "\n")
library(CellChat)

color_data <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/color.csv")
color_mapping <- setNames(color_data$Color, color_data$leiden_mEPC_CAF)

groups_in_network <- levels(cellchat@idents)
color.use <- color_mapping[groups_in_network]
if (any(is.na(color.use))) {
  warning("有些群体标签在 color_mapping 中没有对应的颜色。请检查 color.csv 文件。")
}

png("F2B.png", width = 800, height = 800)

selected_cells <- grep("^(mEPC 08|CAF 06|CAF 04)$", levels(cellchat@idents), value = TRUE)
source_idx <- which(rownames(cellchat@net$weight) %in% selected_cells)
target_idx <- which(colnames(cellchat@net$weight) %in% selected_cells)

mask <- matrix(FALSE, nrow = nrow(cellchat@net$weight), ncol = ncol(cellchat@net$weight))
rownames(mask) <- rownames(cellchat@net$weight)
colnames(mask) <- colnames(cellchat@net$weight)

mask[source_idx, ] <- TRUE
mask[, target_idx] <- TRUE

filtered_net <- cellchat@net$weight
filtered_net[!mask] <- 0

netVisual_circle(
  filtered_net,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  edge.width.max = 15,
  arrow.width = 2,
  arrow.size = 0.5,
  color.use = color.use,
  vertex.label.cex = 0.01
)

dev.off()

png("SF3A.png", width = 800, height = 800)

netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,           
                 weight.scale = TRUE,                 
                 label.edge = FALSE,                  
                 edge.width.max = 15,                 
                 arrow.width = 2,                     
                 arrow.size = 0.5,                    
                 color.use = color.use,               
                 vertex.label.cex = 0.01             
)

dev.off()
## 针对其他两个rds进行类似绘图！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
library(CellChat)

pathways.show <- c("TGFb")
save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF10/figures"

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

png(filename = file.path(save_path, "SF10_TGFB_CellChat.png"), width = 800, height = 800)
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = color.use, vertex.label.cex = 0.001)
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF10/figures"

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

png(filename = file.path(save_path, "SF10_CellChat_TGFB_CellFunction.png"), width = 800, height = 800)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 20,
    height = 8, color.use = color.use, font.size = 14)
dev.off()

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF3/figures"

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = color.use, width = 16, height = 20, font.size = 14)

png(filename = file.path(save_path, "SF3B_outgoing.png"), width = 800, height = 800)
print(ht1)
dev.off()

ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = color.use, width = 16, height = 20, font.size = 14)

png(filename = file.path(save_path, "SF3B_incoming.png"), width = 800, height = 800)
print(ht2)
dev.off()

library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(CellChat)

figures_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B"
if (!dir.exists(figures_path)) {
    dir.create(figures_path, recursive = TRUE)
}
setwd(figures_path)

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/MNFHC_for_cellchat.rds")
head(seurat_obj@meta.data, 5)

seurat_obj@meta.data$subgroup <- ifelse(
  seurat_obj@meta.data$group == "Pr", 
  "Pr", 
  ifelse(
    seurat_obj@meta.data$group == "Me" & seurat_obj@meta.data$batch == "GSE234933", 
    "Dis_Me", 
    ifelse(
      seurat_obj@meta.data$group == "Me" & seurat_obj@meta.data$batch != "GSE234933", 
      "LN_Me", 
      NA
    )
  )
)

seurat_obj@meta.data$group <- seurat_obj@meta.data$subgroup
seurat_obj@meta.data$subgroup <- NULL
head(seurat_obj@meta.data, 5)

data.input <- seurat_obj@assays$RNA@data
meta <- seurat_obj@meta.data

cell.use <- rownames(meta)[meta$group == "Pr" & meta$leiden_mEPC_CAF != "mEPC 10"]
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$leiden_mEPC_CAF)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "leiden_mEPC_CAF")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "leiden_mEPC_CAF")
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

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/cellchat_Pr_withoutmEPC10.rds"
saveRDS(cellchat, file = save_path)
cat("CellChat object has been saved to", save_path, "\n")

cellchat_Pr <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/cellchat_Pr_withoutmEPC10.rds")
cellchat_Me <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2B/cellchat_LN.rds")
object.list <- list(Me = cellchat_Me, Pr = cellchat_Pr)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
ggsave(
  filename = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF10/figures/SF10_cellchatsignal_PrVSMe.png",
  plot = gg1, width = 5, height = 8
)

## nichenet
library(Seurat)
library(tidyverse)
library(dplyr)
options(timeout=600)
library(cowplot)
library(ggplot2)
custom_lib   <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage/nichenetr-2.2.0"
library(nichenetr, lib.loc = custom_lib)
organism = "human"

target_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/reference"
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

files <- c(
  "weighted_networks.rds"                 = "https://zenodo.org/record/3260758/files/weighted_networks.rds",
  "ligand_target_matrix.rds"              = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
  "lr_network.rds"                        = "https://zenodo.org/record/3260758/files/lr_network.rds",
  "weighted_networks_nsga2r_final.rds"    = "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds",
  "ligand_target_matrix_nsga2r_final.rds" = "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
  "lr_network_human_21122021.rds"         = "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
)

for (fname in names(files)) {
  url <- files[fname]
  dest <- file.path(target_dir, fname)
  download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  message("Downloaded ", fname, " to ", dest)
}

lr_network = readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/reference/lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/reference/ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/reference/weighted_networks_nsga2r_final.rds")

seuratObj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
### mEPC 08
Idents(seuratObj)=seuratObj@meta.data[["leiden_mEPC_CAF"]]

nichenet_output = nichenet_seuratobj_aggregate(    
  seurat_obj = seuratObj, 
  receiver = "mEPC 08", 
  condition_colname = "group", condition_oi = "Me", condition_reference = "Pr", 
  sender = c("CAF 04","CAF 06"), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

nichenet_output$ligand_activities
ligand_activities=as.data.frame(nichenet_output$ligand_activities)
head(ligand_activities)

ligand_activities_subset <- ligand_activities %>%
  mutate(test_ligand = reorder(test_ligand, aupr)) %>%
  head(30)

write.csv(
  ligand_activities, 
  file = file.path(out_dir, "mEPC08_ligand_activities_full.csv"), 
  row.names = FALSE
)
write.csv(
  ligand_activities_subset, 
  file = file.path(out_dir, "mEPC08_ligand_activities_top30.csv"), 
  row.names = FALSE
)

min_val <- min(ligand_activities_subset$aupr)
max_val <- max(ligand_activities_subset$aupr)
mid_val <- (min_val + max_val) / 2

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

png(
  filename = file.path(out_dir, "nichenet_ligand_AUPR.png"),
  width    = 3,
  height   = 8,
  units    = "in",
  res      = 300
)

library(ggplot2)
ggplot(ligand_activities_subset, aes(x = factor(1), y = test_ligand, fill = aupr)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "white",
    high   = "orange",
    breaks = c(min_val, mid_val, max_val),
    labels = format(c(min_val, mid_val, max_val), digits = 2)
  ) +
  labs(x = NULL, y = "Prioritized ligands", fill = "AUPR") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

dev.off()
nichenet_output$top_ligands

png(
  filename = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08/nichenet_dotplot.png",
  width    = 8,    
  height   = 6,   
  units    = "in",
  res      = 300   
)
print(nichenet_output$ligand_expression_dotplot)
dev.off()
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
png(
  filename = file.path(output_dir, "ligand_differential_expression_heatmap.png"),
  width    = 1600,      
  height   = 1200,    
  res      = 300       
)
print(nichenet_output$ligand_differential_expression_heatmap)
dev.off()
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
nichenet_output$ligand_target_df %>% head()
nichenet_output$top_targets

output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
TF = as.data.frame(nichenet_output[["ligand_target_heatmap"]][["data"]])
TGFB1 = TF[TF$y == "TGFB1", ]
Ligand = arrange(TGFB1, desc(score))
Top_Ld = as.character(Ligand$x[1:20])
top_targets <- TF %>%
  filter(y == "TGFB1") %>%
  arrange(desc(score)) %>%
  slice_head(n = 20) %>%
  pull(x)

top_targets <- as.character(top_targets)
output_file <- file.path(output_dir, "Top_ligand_targets_TGFB1.txt")
writeLines(top_targets, output_file)
message("已将前20个 ligand 写入：", output_file)
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_path <- file.path(output_dir, "nichenet_heatmap_mepc08_polished.png")
png(
  filename = output_path,
  width    = 12,
  height   = 8,
  units    = "in",
  res      = 300
)
print(
  nichenet_output$ligand_target_heatmap +
    scale_fill_gradient2(low  = "whitesmoke", high = "royalblue") +
    xlab("Predicted target genes in Malignant mEPC 08") +
    ylab("Prioritized CAF 04 and CAF 06 ligands") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
    )
)
dev.off()
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08"
output_file <- file.path(output_dir, "nichenet_target_dotplot_mepc08.png")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

mal8_subset <- subset(seuratObj, idents = "mEPC 08")
p <- DotPlot(
  object    = mal8_subset,
  features  = Top_Ld,
  cols      = "RdYlBu",
  split.by  = "group"
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(
  filename = output_file,
  plot     = p,
  width    = 10,
  height   = 5,
  units    = "in",
  dpi      = 300
)

Mal8_Seurat = seuratObj %>% subset(idents = "mEPC 08")
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08/nichenet_target_vlnplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
mal8_subset <- subset(seuratObj, idents = "mEPC 08")

output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08/nichenet_target_vlnplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
mal8_subset <- subset(seuratObj, idents = "mEPC 08")
group_colors <- c("#F8766D", "#078992")
for (gene in Top_Ld) {
  p <- VlnPlot(
    object   = mal8_subset,
    features = gene,
    split.by = "group",
    pt.size  = 0,
    cols     = group_colors
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(gene)
  
  out_file <- file.path(output_dir, paste0(gene, "_vlnplot.png"))
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 4,
    height   = 4,
    units    = "in",
    dpi      = 300
  )
}

output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/mEPC08/nichenet_target_vlnplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
mal8_subset <- subset(seuratObj, idents = "mEPC 08")
group_colors <- c("#F8766D", "#078992")  
genes <- c("TGFBR1", "TGFBR2", "TGFBR3", "WDR54")

# 5. 循环绘图并保存
for (gene in genes) {
  p <- VlnPlot(
    object   = mal8_subset,
    features = gene,
    split.by = "group",
    pt.size  = 0,
    cols     = group_colors
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(gene)
  
  out_file <- file.path(output_dir, paste0(gene, "_vlnplot.png"))
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 4,
    height   = 4,
    units    = "in",
    dpi      = 300
  )
}

## CAF 06
seuratObj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
Idents(seuratObj)=seuratObj@meta.data[["leiden_mEPC_CAF"]]
nichenet_output2 = nichenet_seuratobj_aggregate(
  seurat_obj             = seuratObj, 
  receiver               = "CAF 06", 
  condition_colname      = "group", 
  condition_oi           = "Me", 
  condition_reference    = "Pr", 
  sender                 = c("CAF 04","mEPC 08"), 
  ligand_target_matrix   = ligand_target_matrix,
  lr_network             = lr_network,
  weighted_networks      = weighted_networks
)
nichenet_output2$ligand_activities
ligand_activities = as.data.frame(nichenet_output2$ligand_activities)
head(ligand_activities)
ligand_activities_subset <- ligand_activities %>%
  mutate(test_ligand = reorder(test_ligand, aupr)) %>%
  head(30)


write.csv(
  ligand_activities, 
  file = file.path(out_dir, "CAF06_ligand_activities_full.csv"), 
  row.names = FALSE
)

write.csv(
  ligand_activities_subset, 
  file = file.path(out_dir, "CAF06_ligand_activities_top30.csv"), 
  row.names = FALSE
)

min_val <- min(ligand_activities_subset$aupr)
max_val <- max(ligand_activities_subset$aupr)
mid_val <- (min_val + max_val) / 2
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/CAF06"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
p <- ggplot(ligand_activities_subset, aes(
      x    = factor(1),
      y    = test_ligand,
      fill = aupr
    )) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "white",
    high   = "orange",
    breaks = c(min_val, mid_val, max_val),
    labels = format(c(min_val, mid_val, max_val), digits = 2)
  ) +
  labs(
    x    = "",
    y    = "Prioritized ligands",
    fill = "AUPR"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )
out_file <- file.path(output_dir, "nichenet_ligand_AUPR_CAF06.png")
ggsave(
  filename = out_file,
  plot     = p,
  width    = 3,    
  height   = 8,   
  units    = "in",
  dpi      = 300
)
nichenet_output2$top_ligands
png(
  filename = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/CAF06/nichenet_dotplot.png",
  width    = 8,   
  height   = 6,    
  units    = "in",
  res      = 300   
)
print(nichenet_output2$ligand_expression_dotplot)
dev.off()
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/CAF06"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
out_file <- file.path(output_dir, "nichenet_ligand_DE_heatmap.png")
png(
  filename = out_file,
  width    = 9,    
  height   = 8,    
  units    = "in",
  res      = 300
)
print(nichenet_output2$ligand_differential_expression_heatmap)
dev.off()
nichenet_output2$ligand_target_matrix %>% .[1:10, 1:6]
nichenet_output2$ligand_target_df %>% head()
nichenet_output2$top_targets
TF      = as.data.frame(nichenet_output2[["ligand_target_heatmap"]][["data"]])
TGFB1   = TF[TF$y == "TGFB1", ]
Ligand  = arrange(TGFB1, desc(score))
Top_Ld  = Ligand$x[1:20]
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/SF9/nichnet/figure/CAF06"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
caf06_subset  <- subset(seuratObj, idents = "CAF 06")
group_colors  <- c("#F8766D", "#078992")  # Me = rgb(7,137,146), Pr = rgb(248,118,109)
png(
  filename = file.path(output_dir, "ligand_target_heatmap_CAF06.png"),
  width    = 12, height = 8, units = "in", res = 300
)
print(
  nichenet_output2$ligand_target_heatmap +
    scale_fill_gradient2(low = "whitesmoke", high = "royalblue") +
    xlab("Predicted target genes in CAF 06 ligands") +
    ylab("Prioritized CAF 04 and Malignant mEPC 08")
)
dev.off()
png(
  filename = file.path(output_dir, "ligand_DE_heatmap_CAF06.png"),
  width    = 9, height = 8, units = "in", res = 300
)
print(nichenet_output2$ligand_differential_expression_heatmap)
dev.off()
png(
  filename = file.path(output_dir, "ligand_expression_dotplot_CAF06.png"),
  width    = 10, height = 5, units = "in", res = 300
)
print(nichenet_output2$ligand_expression_dotplot)
dev.off()
png(
  filename = file.path(output_dir, "ligand_activity_target_heatmap_CAF06.png"),
  width    = 12, height = 8, units = "in", res = 300
)
print(nichenet_output2$ligand_activity_target_heatmap)
dev.off()
png(
  filename = file.path(output_dir, "ligand_receptor_heatmap_CAF06.png"),
  width    = 9, height = 8, units = "in", res = 300
)
print(nichenet_output2$ligand_receptor_heatmap)
dev.off()
png(
  filename = file.path(output_dir, "top_ligands_dotplot_CAF06.png"),
  width    = 10, height = 5, units = "in", res = 300
)
print(
  DotPlot(
    caf06_subset,
    features = Top_Ld,
    cols     = "RdYlBu",
    split.by = "group"
  ) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

for (gene in Top_Ld) {
  png(
    filename = file.path(output_dir, paste0(gene, "_vlnplot_CAF06.png")),
    width    = 4, height = 4, units = "in", res = 300
  )
  print(
    VlnPlot(
      caf06_subset,
      features = gene,
      split.by = "group",
      pt.size  = 0,
      cols     = group_colors
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(gene)
  )
  dev.off()
}
