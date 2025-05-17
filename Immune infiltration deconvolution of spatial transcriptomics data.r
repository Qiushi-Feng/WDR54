#### ST-IOBR
#### This script serves only as an example for a single sample; in practice, you should apply the same code to the other samples.
library(tidyverse)
library(IOBR)
library(reshape2)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(cowplot)
library(ComplexHeatmap)
library(SPOTlight)
library(SingleCellExperiment)
library(magick)
library(SpatialExperiment)
library(scater)
library(scran)
library(scatterpie)
library(patchwork)
library(Seurat)
library(rtracklayer)
library(GenomicRanges)
library(Matrix)
library(SpatialExperiment)

gtf_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/Rpackage/Homo_sapiens.GRCh38.113.chr.gtf"
cat("[INFO] 读取 GTF 并计算基因长度 …\n")
gtf <- import(gtf_file)
exon_gr <- gtf[gtf$type == "exon"]
grl <- split(exon_gr, exon_gr$gene_id)
gene_len_bp <- vapply(grl, function(gr) sum(width(reduce(gr))), numeric(1))
gene_symbol_map <- unique(
  data.frame(
    gene_id   = exon_gr$gene_id,
    gene_name = exon_gr$gene_name,
    stringsAsFactors = FALSE
  )
)
symbol2id <- setNames(gene_symbol_map$gene_id, gene_symbol_map$gene_name)

in_dir  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata"
out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (i in 1:4) {
  sample_id <- paste0("Patient", i)
  rds_path  <- file.path(in_dir, paste0(sample_id, "_processed.rds"))
  out_tsv   <- file.path(out_dir, paste0(sample_id, "_SCT_TPM.tsv"))
  
  if (!file.exists(rds_path)) {
    warning("[WARN] 文件不存在：", rds_path)
    next
  }
  cat("[INFO] 处理 ", sample_id, " …\n")
  
  obj <- readRDS(rds_path)
  if (!"SCT" %in% Assays(obj)) {
    warning(sample_id, " 缺少 SCT assay，跳过。")
    next
  }
  count_mat <- GetAssayData(obj[["SCT"]], slot = "counts")
  
  genes   <- rownames(count_mat)
  len_vec <- gene_len_bp[genes]
  
  need_sym <- is.na(len_vec)
  if (any(need_sym)) {
    matched_id        <- symbol2id[genes[need_sym]]
    len_vec[need_sym] <- gene_len_bp[matched_id]
  }
  keep <- !is.na(len_vec)
  if (sum(!keep) > 0) {
    warning(sum(!keep), " 个基因未匹配到长度，已剔除。")
  }
  count_mat <- count_mat[keep, ]
  len_kb    <- len_vec[keep] / 1e3
  
  rpk         <- sweep(count_mat, 1, len_kb, FUN = "/")
  col_sum_rpk <- Matrix::colSums(rpk)
  tpm         <- sweep(rpk, 2, col_sum_rpk / 1e6, FUN = "/")
  
  write.table(
    as.matrix(tpm),
    file         = out_tsv,
    sep          = "\t",
    quote        = FALSE,
    row.names    = TRUE,
    col.names    = NA,
    fileEncoding = "UTF-8"
  )
  cat("    ↳ 已写出：", out_tsv, "\n")
}

cat("[INFO] 全部样本处理完毕（保留全部基因，输出 TSV）。\n")

input_dir   <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/data"
result_root <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/result"

for (i in 1:4) {
  patient    <- paste0("Patient", i)
  data_file  <- file.path(input_dir, paste0(patient, "_SCT_TPM.tsv"))
  result_dir <- file.path(result_root, patient)

  data <- read_tsv(
    data_file,
    col_names      = TRUE,
    show_col_types = FALSE
  ) %>%
    column_to_rownames(var = "...1") %>%
    as.data.frame()

  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  cibersort <- deconvo_tme(
    eset   = data,
    method = "cibersort",
    arrays = FALSE,
    perm   = 50
  )
  write_csv(cibersort, file.path(result_dir, "cibersort.csv"))

  epic <- deconvo_tme(
    eset   = data,
    method = "epic",
    arrays = FALSE
  )
  write_csv(epic, file.path(result_dir, "epic.csv"))

  mcp <- deconvo_tme(
    eset   = data,
    method = "mcpcounter"
  )
  write_csv(mcp, file.path(result_dir, "mcpcounter.csv"))

  xcell <- deconvo_tme(
    eset   = data,
    method = "xcell",
    arrays = FALSE
  )
  write_csv(xcell, file.path(result_dir, "xcell.csv"))

  timer <- deconvo_tme(
    eset       = data,
    method     = "timer",
    group_list = rep("uvm", ncol(data))
  )
  write_csv(timer, file.path(result_dir, "timer.csv"))

  quantiseq <- deconvo_tme(
    eset       = data,
    method     = "quantiseq",
    tumor      = TRUE,
    arrays     = FALSE,
    scale_mrna = TRUE
  )
  write_csv(quantiseq, file.path(result_dir, "quantiseq.csv"))

  estimate_result <- deconvo_tme(
    eset   = data,
    method = "estimate"
  )
  write_csv(estimate_result, file.path(result_dir, "estimate.csv"))

  ips_result <- deconvo_tme(
    eset   = data,
    method = "ips",
    plot   = FALSE
  )
  write_csv(ips_result, file.path(result_dir, "ips.csv"))
}

rds_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata/Patient1_processed.rds"
seu <- readRDS(rds_path)

new_barcodes <- sub("^Patient1_", "", colnames(seu))
colnames(seu) <- new_barcodes

sce <- as.SingleCellExperiment(seu, assay = NULL)

coord_csv <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/tissue_positions_list.csv"
coord_df  <- read.csv(coord_csv, header = FALSE, stringsAsFactors = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "row", "col", "imagecol", "imagerow")

coord_df <- coord_df[match(colnames(sce), coord_df$barcode), ]
if (any(is.na(coord_df$barcode))) stop("部分 barcodes 在坐标文件里找不到，请检查前缀/文件路径！")

xy <- as.matrix(coord_df[, c("imagecol", "imagerow")])
rownames(xy) <- coord_df$barcode
colnames(xy) <- c("x", "y")

sf_path       <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/scalefactors_json.json"
sf            <- jsonlite::read_json(sf_path)
lowres_scale  <- sf$tissue_lowres_scalef

png_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/tissue_hires_image.png"
png_img  <- png::readPNG(png_path)
sp_img   <- SpatialImage(x = as.raster(png_img))

imgData <- DataFrame(
  sample_id   = "Patient1",
  image_id    = "slice1",
  data        = I(list(sp_img)),
  scaleFactor = lowres_scale
)

spe <- SpatialExperiment(
  assays        = assays(sce),
  rowData       = rowData(sce),
  colData       = colData(sce),
  spatialCoords = xy,
  sample_id     = rep("Patient1", ncol(sce)),
  imgData       = imgData,
  reducedDims   = reducedDims(sce),
  altExps       = altExps(sce),
  metadata      = metadata(sce)
)

validObject(spe)
coords <- spatialCoords(spe)
coords <- coords * 3.33
spatialCoords(spe) <- coords

prop_path   <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/result/Patient1/quantiseq.csv"
output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/figure/Patient1_quantiseq_pie.png"

mat <- read_csv(prop_path, show_col_types = FALSE) %>%
  column_to_rownames(var = "ID") %>%
  rename_with(~ str_remove(., "_quantiseq$"))

head(mat)
mat[mat < 0.001] <- 0
rownames(mat) <- sub("^Patient1_", "", rownames(mat))
cell_types <- colnames(mat)

custom_colors <- c(
  "#E5AF44", "#44AFE5", "#88CEE6", "#44E5AF", "#B83945",
  "#4D4D4D", "#BC8F8F", "#E89DA0", "#B696B6", "#FB8C62", "#D3D3D3"
)

if (length(cell_types) > length(custom_colors)) {
  warning("细胞类型数量超过颜色数量，颜色将被重复使用。")
}
pal <- rep(custom_colors, length.out = length(cell_types))
names(pal) <- cell_types
p <- plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = cell_types,
  img = TRUE,
  slice = "Patient1",
  scatterpie_alpha = 1,
  pie_scale = 0.4
) +
  scale_fill_manual(values = pal, breaks = names(pal)) +
  ggtitle("Patient1 – quanTIseq") +
  theme_minimal()

ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300)

rds_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/processeddata/Patient1_processed.rds"
seu <- readRDS(rds_path)

new_barcodes <- sub("^Patient1_", "", colnames(seu))
colnames(seu) <- new_barcodes

sce <- as.SingleCellExperiment(seu, assay = NULL)

coord_csv <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/tissue_positions_list.csv"
coord_df <- read.csv(coord_csv, header = FALSE, stringsAsFactors = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "row", "col", "imagecol", "imagerow")

coord_df <- coord_df[match(colnames(sce), coord_df$barcode), ]
if (any(is.na(coord_df$barcode))) stop("部分 barcodes 在坐标文件里找不到，请检查前缀/文件路径！")

xy <- as.matrix(coord_df[, c("imagecol", "imagerow")])
rownames(xy) <- coord_df$barcode
colnames(xy) <- c("x", "y")

sf_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/scalefactors_json.json"
sf <- jsonlite::read_json(sf_path)
lowres_scale <- sf$tissue_lowres_scalef

png_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2_6/rawdata/GSE220978_RAW/image/Patient1_image/tissue_hires_image.png"
png_img <- png::readPNG(png_path)
sp_img <- SpatialImage(x = as.raster(png_img))

imgData <- DataFrame(
  sample_id = "Patient1",
  image_id = "slice1",
  data = I(list(sp_img)),
  scaleFactor = lowres_scale
)

spe <- SpatialExperiment(
  assays = assays(sce),
  rowData = rowData(sce),
  colData = colData(sce),
  spatialCoords = xy,
  sample_id = rep("Patient1", ncol(sce)),
  imgData = imgData,
  reducedDims = reducedDims(sce),
  altExps = altExps(sce),
  metadata = metadata(sce)
)

validObject(spe)
coords <- spatialCoords(spe)
coords <- coords * 3.33
spatialCoords(spe) <- coords

prop_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/result/Patient1/cibersort.csv"
output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/figure/Patient1_CIBERSORT_pie.png"

mat <- read_csv(prop_path, show_col_types = FALSE) %>%
  column_to_rownames(var = "ID") %>%
  rename_with(~ str_remove(., "_CIBERSORT$"))

mat <- mat[, -((ncol(mat) - 2):ncol(mat))]

mat <- sweep(mat, 1, rowSums(mat), FUN = "/")

head(mat)
mat[mat < 0.001] <- 0
rownames(mat) <- sub("^Patient1_", "", rownames(mat))
cell_types <- colnames(mat)

custom_colors <- c(
  "#A25A59", "#993F47", "#66484B", "#4CAFDD", "#61BCE5",
  "#82CEE6", "#4C9BE6", "#6ED7D1", "#4EE2B7", "#6BAC8C",
  "#675D5D", "#9C7C7C", "#D9989A", "#E19CA3", "#C999AE",
  "#BD95AE", "#E5AF44", "#F98F67", "#4D4D4D", "#D3D3D3",
  "#EF607A", "#6BBC47"
)

if (length(cell_types) > length(custom_colors)) {
  warning("细胞类型数量超过颜色数量，颜色将被重复使用。")
}
pal <- rep(custom_colors, length.out = length(cell_types))
names(pal) <- cell_types

p <- plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = cell_types,
  img = FALSE,
  slice = "Patient1",
  scatterpie_alpha = 1,
  pie_scale = 0.4
) +
  scale_fill_manual(values = pal, breaks = names(pal)) +
  ggtitle("Patient1 – CIBERSORT") +
  theme_minimal()

ggsave(output_path, plot = p, width = 15, height = 12, dpi = 300)

input_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F4/ST-IOBR/figure/Patient1_CIBERSORT_pie.png"

img <- image_read(input_path)
img <- image_rotate(img, 180)
img <- image_flop(img)

image_write(img, path = input_path)

