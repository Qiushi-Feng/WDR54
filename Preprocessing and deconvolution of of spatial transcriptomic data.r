### Object construction and preprocessing
library(Seurat)
library(Matrix)
library(dplyr)
library(harmony)
library(RColorBrewer)
library(ggplot2)
library(SpatialExperiment)
library(monocle)
library(Seurat)
library(argparse)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
library(scater)
library(data.table)
library(scran)
library(SPOTlight)
library(SingleCellExperiment)
library(tidyverse)
library(SeuratDisk)

patient_data_dir  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/rawdata/GSE220978_RAW/data/Patient1"
patient_image_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/rawdata/GSE220978_RAW/image/Patient1_image"
matrix_file   <- file.path(patient_data_dir, "GSM6833484_Patient1_matrix.mtx")
barcodes_file <- file.path(patient_data_dir, "GSM6833484_Patient1_barcodes.tsv")
features_file <- file.path(patient_data_dir, "GSM6833484_Patient1_features.tsv")
mat <- readMM(matrix_file)
barcodes <- readLines(barcodes_file)
features_data <- read.table(features_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene_names <- if (ncol(features_data) >= 2) features_data[, 2] else features_data[, 1]
rownames(mat) <- gene_names
colnames(mat) <- barcodes
seurat_obj <- CreateSeuratObject(
  counts = mat,
  project = "Patient1",
  assay = "Spatial"
)
img <- Read10X_Image(
  image.dir  = patient_image_dir,
  image.name = "tissue_hires_image.png",
  filter.matrix = TRUE
)
DefaultAssay(img) <- "Spatial"
seurat_obj[["slice1"]] <- img
rds_file <- file.path(patient_data_dir, "Patient1_spatial.rds")
saveRDS(seurat_obj, file = rds_file)
patients <- c("Patient1", "Patient2", "Patient3", "Patient4")
raw_base_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/rawdata/GSE220978_RAW/data"
processed_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
  cat("[INFO] 创建目录:", processed_dir, "\n")
}
for (patient in patients) {
  file_path <- file.path(raw_base_dir, patient, paste0(patient, "_spatial.rds"))
  sample_object <- readRDS(file_path)
  cat("[INFO] 已读取文件:", file_path, "\n")
  sample_object@project.name <- patient
  sample_object@meta.data$orig.ident <- patient
  sample_object <- RenameCells(object = sample_object, add.cell.id = patient)
  sample_object <- SCTransform(sample_object, assay = "Spatial", method = "poisson")
  cat("[INFO] SCTransform 已完成: ", patient, "\n")
  output_file <- file.path(processed_dir, paste0(patient, "_SCT.rds"))
  saveRDS(sample_object, file = output_file)
  cat("[INFO] 保存后的文件已输出至:", output_file, "\n")
}
patients <- c("Patient1", "Patient2", "Patient3", "Patient4")
input_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
figure_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures"
output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("[INFO] 创建输出目录:", output_dir, "\n")
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
  cat("[INFO] 创建图像目录:", figure_dir, "\n")
}
scale_factor <- 3.33
slice_names <- c("slice1", "slice2", "slice3", "slice4")
for (patient in patients) {
  file_path <- file.path(input_dir, paste0(patient, "_SCT.rds"))
  ST <- readRDS(file_path)
  cat("[INFO]", patient, "ST 对象已成功读取：", file_path, "\n")
  ST <- RunPCA(ST, assay = "SCT", features = VariableFeatures(ST), verbose = FALSE)
  cat("[INFO]", patient, "PCA 降维完成。\n")
  ST <- FindNeighbors(ST, reduction = "pca", dims = 1:30, k.param = 30) %>%
        FindClusters(resolution = 1.5, algorithm = 4, group.singletons = TRUE)
  cat("[INFO]", patient, "邻域构建与聚类完成。\n")
  cat("[INFO]", patient, "每个聚类中不同 orig.ident 的细胞数量：\n")
  cluster_origin_counts <- table(ST@meta.data$seurat_clusters, ST@meta.data$orig.ident)
  print(cluster_origin_counts)
  ST <- RunUMAP(ST, reduction = "pca", dims = 1:30, n.neighbors = 30)
  cat("[INFO]", patient, "UMAP 降维完成。\n")
  plot_file <- file.path(figure_dir, paste0(patient, "_umap.png"))
  png(filename = plot_file, width = 2400, height = 1800, res = 300)
  p <- DimPlot(ST, reduction = "umap", pt.size = 2)
  p$layers[[1]]$aes_params$alpha <- 0.5  
  p$layers[[1]]$aes_params$shape <- 16
  print(p)
  dev.off()
  cat("[INFO]", patient, "UMAP 图像已保存至:", plot_file, "\n\n")
  for (slice in slice_names) {
    if (!is.null(ST@images[[slice]])) {
      img_obj <- ST@images[[slice]]
      coords <- img_obj@coordinates
      coords$imagecol <- coords$imagecol * scale_factor
      coords$imagerow <- coords$imagerow * scale_factor
      img_obj@coordinates <- coords
      ST@images[[slice]] <- img_obj
      cat(sprintf("[INFO] 图像 %s 坐标已放大 %.2f 倍。\n", slice, scale_factor))
    } else {
      cat(sprintf("[INFO] 图像 %s 不存在，跳过。\n", slice))
    }
  }
  output_file <- file.path(output_dir, paste0(patient, "_processed.rds"))
  saveRDS(ST, file = output_file)
  cat("[INFO]", patient, "处理后的 ST 对象已保存至:", output_file, "\n\n")
}
library(Seurat)
in_dir  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/figures"
pal <- list(
  Patient1 = c(`1`="#2171A8",`2`="#EC7C1E",`3`="#37955F",`4`="#642D84",`5`="#D74692",
               `6`="#EADE25",`7`="#8C564B",`8`="#6AC3CA",`9`="#C5B0D5",`10`="#F4B579",
               `11`="#98CC85",`12`="#D62728",`13`="#B2DF8A",`14`="#BABABA"),
  Patient2 = c(`1`="#2171A8",`2`="#EC7C1E",`3`="#37955F",`4`="#642D84",`5`="#C5B0D5",
               `6`="#D74692",`7`="#EADE25",`8`="#D62728",`9`="#6AC3CA",`10`="#F4B579",
               `11`="#8C564B",`12`="#98CC85",`13`="#B2DF8A",`14`="#BABABA",
               `15`="#2F2D54",`16`="#E0A4DD"),
  Patient3 = c(`1`="#2171A8",`2`="#EC7C1E",`3`="#D62728",`4`="#C5B0D5",`5`="#8C564B",
               `6`="#37955F",`7`="#642D84",`8`="#D74692",`9`="#EADE25",`10`="#6AC3CA",
               `11`="#F4B579",`12`="#98CC85"),
  Patient4 = c(`1`="#2171A8",`2`="#C5B0D5",`3`="#EC7C1E",`4`="#37955F",`5`="#642D84",
               `6`="#D62728",`7`="#8C564B",`8`="#D74692",`9`="#EADE25",`10`="#6AC3CA",
               `11`="#F4B579",`12`="#98CC85",`13`="#B2DF8A")
)
for(i in 1:4){
  pid   <- paste0("Patient", i)
  obj   <- readRDS(file.path(in_dir, paste0(pid, "_processed.rds")))
  idlev <- levels(factor(obj$seurat_clusters))
  cols  <- pal[[pid]][idlev]
  cols[is.na(cols)] <- "#CCCCCC"
  png(file.path(out_dir, paste0(pid, "_cluster.png")), width = 3000, height = 3000, res = 300)
  print(SpatialPlot(obj, group.by = "seurat_clusters", cols = cols, crop = FALSE, pt.size.factor = 6.5, stroke = 0, label = TRUE, label.size = 7, label.color = "black", alpha = 0.8) + theme(legend.position = "none"))
  dev.off()
  cat("Saved:", pid, "\n")
}

### Deconvolution
obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
obj <- SCTransform(obj, assay = "RNA")
DefaultAssay(obj) <- "SCT"
cat("[INFO] SCTransform 处理完成，SCT assay 已设为默认\n")
saveRDS(obj, "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds")
scRNA <- readRDS('/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2D/vector/MNFHC_for_vector.rds')
sce <- as.SingleCellExperiment(scRNA)
idx <- split(seq(ncol(sce)), sce$leiden_mEPC_CAF)
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if(n < 150) {
    sample_count <- n
  } else {
    sample_count <- max(round(n * 0.1), 150)
  }
  cat("组:", unique(sce$leiden_mEPC_CAF[i]),
      "原始细胞数:", n,
      "抽取后细胞数:", sample_count, "\n")
  sample(i, sample_count)
})
sce <- sce[, unlist(cs_keep)]
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sce) <- colData(sce)$leiden_mEPC_CAF
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
mgs <- scoreMarkers(sce, subset.row = genes)
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  x <- x[x$mean.AUC > 0.7, ]
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient1_processed.rds")
stRNA <- as.SingleCellExperiment(seurat_obj)
res <- SPOTlight(
  x = sce,
  y = stRNA,
  groups = sce$leiden_mEPC_CAF,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
save(res, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient1_processed.rda")
head(mat <- res$mat)
mod <- res$NMF
res <- get(load("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient1_processed.rda"))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient1_processed.rds")
head(mat <- res$mat)
mod <- res$NMF
seurat_obj[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(seurat_obj) <- "SPOTlight"
features <- rownames(seurat_obj[["SPOTlight"]])
for(feature in features) {
  output_file <- paste0(
    "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/projection/Patient1_",
    feature, ".png"
  )
  png(filename = output_file, width = 3000, height = 3000, res = 300)
  p <- SpatialFeaturePlot(seurat_obj, 
                          features = feature, 
                          pt.size.factor = 6.5, 
                          stroke = 0, 
                          alpha = 1, 
                          crop = FALSE)
  print(p)
  dev.off()
}

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient2_processed.rds")
stRNA <- as.SingleCellExperiment(seurat_obj)
res <- SPOTlight(
  x = sce,
  y = stRNA,
  groups = sce$leiden_mEPC_CAF,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
save(res, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient2_processed.rda")
head(mat <- res$mat)
mod <- res$NMF
res <- get(load("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient2_processed.rda"))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient2_processed.rds")
head(mat <- res$mat)
mod <- res$NMF
seurat_obj[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(seurat_obj) <- "SPOTlight"
features <- rownames(seurat_obj[["SPOTlight"]])
for(feature in features) {
  output_file <- paste0(
    "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/projection/Patient2_",
    feature, ".png"
  )
  png(filename = output_file, width = 3000, height = 3000, res = 300)
  p <- SpatialFeaturePlot(seurat_obj, 
                          features = feature, 
                          pt.size.factor = 6.5, 
                          stroke = 0, 
                          alpha = 1, 
                          crop = FALSE)
  print(p)
  dev.off()
}

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient3_processed.rds")
stRNA <- as.SingleCellExperiment(seurat_obj)
res <- SPOTlight(
  x = sce,
  y = stRNA,
  groups = sce$leiden_mEPC_CAF,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
save(res, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient3_processed.rda")
head(mat <- res$mat)
mod <- res$NMF
res <- get(load("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient3_processed.rda"))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient3_processed.rds")
head(mat <- res$mat)
mod <- res$NMF
seurat_obj[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(seurat_obj) <- "SPOTlight"
features <- rownames(seurat_obj[["SPOTlight"]])
for(feature in features) {
  output_file <- paste0(
    "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/projection/Patient3_",
    feature, ".png"
  )
  png(filename = output_file, width = 3000, height = 3000, res = 300)
  p <- SpatialFeaturePlot(seurat_obj, 
                          features = feature, 
                          pt.size.factor = 6.5, 
                          stroke = 0, 
                          alpha = 1, 
                          crop = FALSE)
  print(p)
  dev.off()
}

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient4_processed.rds")
stRNA <- as.SingleCellExperiment(seurat_obj)
res <- SPOTlight(
  x = sce,
  y = stRNA,
  groups = sce$leiden_mEPC_CAF,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
save(res, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient4_processed.rda")
head(mat <- res$mat)
mod <- res$NMF
res <- get(load("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/Patient4_processed.rda"))
seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/Patient4_processed.rds")
head(mat <- res$mat)
mod <- res$NMF
seurat_obj[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(seurat_obj) <- "SPOTlight"
features <- rownames(seurat_obj[["SPOTlight"]])
for(feature in features) {
  output_file <- paste0(
    "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight/projection/Patient4_",
    feature, ".png"
  )
  png(filename = output_file, width = 3000, height = 3000, res = 300)
  p <- SpatialFeaturePlot(seurat_obj, 
                          features = feature, 
                          pt.size.factor = 6.5, 
                          stroke = 0, 
                          alpha = 1, 
                          crop = FALSE)
  print(p)
  dev.off()
}
base_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight"
files <- paste0(base_path, "/Patient", 1:4, "_processed.rda")
for (f in files) {
  load(f)
  df <- as.data.frame(res$mat, stringsAsFactors = FALSE)
  df$ID <- rownames(res$mat)
  sample_name <- sub("_processed\\.rda", "", basename(f))
  df$Sample <- sample_name
  output_file <- paste0(base_path, "/", sample_name, "_res_mat.csv")
  fwrite(df, file = output_file)
  cat("样本", sample_name, "的数据已保存到:", output_file, "\n")
}
base_csv_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight"
base_rda_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata"
for(i in 1:4) {
  csv_file <- paste0(base_csv_path, "/Patient", i, "_res_mat.csv")
  rda_file <- paste0(base_rda_path, "/Patient", i, "_processed.rds")
  csv_data <- read.csv(csv_file, stringsAsFactors = FALSE)
  loaded_objs <- load(rda_file)
  st_obj <- get(loaded_objs[1])
  cell_order <- rownames(st_obj@meta.data)
  csv_data_sorted <- csv_data[match(cell_order, csv_data$ID), ]
  csv_data_sorted$seurat_clusters <- st_obj@meta.data[cell_order, "seurat_clusters"]
  write.csv(csv_data_sorted, file = csv_file, row.names = FALSE)
  cat("Patient", i, ": CSV文件已根据 RDA 中细胞索引顺序排序，并成功添加了 seurat_clusters 列。\n")
}
base_csv_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spotlight/SCTobjectivespotlight"
final_result <- list()
for(i in 1:4) {
  csv_file <- file.path(base_csv_path, paste0("Patient", i, "_res_mat.csv"))
  csv_data <- read.csv(csv_file, stringsAsFactors = FALSE)
  features <- grep("^(CAF|mEPC)", names(csv_data), value = TRUE)
  result_list <- list()
  for(feat in features) {
    avg_vals <- tapply(csv_data[[feat]], csv_data$seurat_clusters, mean, na.rm = TRUE)
    best_cluster <- names(which.max(avg_vals))
    best_mean <- max(avg_vals, na.rm = TRUE)
    result_list[[feat]] <- data.frame(
      Patient = paste0("Patient", i),
      Feature = feat,
      BestCluster = best_cluster,
      MeanValue = best_mean,
      stringsAsFactors = FALSE
    )
  }
  result_df <- do.call(rbind, result_list)
  final_result[[paste0("Patient", i)]] <- result_df
  cat("---- Results for Patient", i, "----\n")
  print(result_df)
}
all_results <- do.call(rbind, final_result)
output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/processeddata/spolight_results.csv"
write.csv(all_results, file = output_file, row.names = FALSE)
cat("==== Combined Results have been saved to:", output_file, "\n")
