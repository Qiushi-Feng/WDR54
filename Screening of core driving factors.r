### Screening of core driving factors
##  hdWGCNA(mEPC 08)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(UCell)

work_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3A"
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
  cat("Directory created:", work_dir, "\n")
} else {
  cat("Directory already exists:", work_dir, "\n")
}
setwd(work_dir)
cat("Working directory set to:", getwd(), "\n")

theme_set(theme_cowplot())
set.seed(123)

seurat_obj <- readRDS("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3A/MNFHC_mEPC.rds")

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial"
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("leiden_mEPC_CAF"),
  reduction = 'harmony',
  ident.group = 'leiden_mEPC_CAF'
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "mEPC 08",
  group.by = 'leiden_mEPC_CAF',
  assay = 'RNA',
  slot = 'data'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

plot_list <- PlotSoftPowers(seurat_obj)

output_file <- "SF4A.png"
png(filename = output_file, width = 1200, height = 800, res = 150)
wrap_plots(plot_list, ncol = 2)
dev.off()

cat("软阈值评估图已保存到", output_file, "\n")

power_table <- GetPowerTable(seurat_obj)
head(power_table)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'mEPC 08'
)

output_file <- "F3A_1.png"
png(filename = output_file, width = 1200, height = 800, res = 150)
PlotDendrogram(seurat_obj, main = NULL)
dev.off()

cat("Dendrogram已保存为:", output_file, "\n")

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "leiden_mEPC_CAF"
)

hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized = FALSE)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'leiden_mEPC_CAF',
  group_name = 'mEPC 08'
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "mEPC 08-M"
)

p <- PlotKMEs(seurat_obj, ncol = 7)

ggsave(filename = "SF4B_1.png", plot = p, path = getwd(), device = "png", width = 16, height = 6)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
head(modules[, 1:10])

save_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3A/modules.csv"
write.csv(modules, file = save_path, row.names = FALSE)
message("Modules 数据已保存至: ", save_path)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 'all',
  method = 'UCell'
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs',
  order = TRUE
)

combined_plot <- wrap_plots(plot_list, ncol = 7)
ggsave("SF4B_2.png", plot = combined_plot, width = 20, height = 10, dpi = 300)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'scores',
  order = 'shuffle',
  ucell = TRUE
)

combined_plot <- wrap_plots(plot_list, ncol = 7)
ggsave("SF4B_3.png", plot = combined_plot, width = 20, height = 10, dpi = 300)

invisible({
  png("SF4C.png", width = 8, height = 6.4, units = "in", res = 300)
  ModuleCorrelogram(seurat_obj)
  dev.off()
})

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(seurat_obj, features = mods, group.by = 'leiden_mEPC_CAF')
p <- p + RotatedAxis() + scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue')
ggsave("SF4B_4.png", plot = p, width = 10, height = 8, dpi = 300)

color_palette <- c(
  'mEPC 01' = '#1F77B4',
  'mEPC 02' = '#FF7F0E',
  'mEPC 03' = '#2CA02C',
  'mEPC 04' = '#D62728',
  'mEPC 05' = '#9467BD',
  'mEPC 06' = '#8C564B',
  'mEPC 07' = '#E377C2',
  'mEPC 08' = '#7F7F7F',
  'mEPC 09' = '#BCBD22',
  'mEPC 10' = '#17BECF'
)

p <- VlnPlot(
  seurat_obj,
  features = 'mEPC 08-M5',
  group.by = 'leiden_mEPC_CAF',
  pt.size = 0
)
p <- p + geom_boxplot(width = 0.25, fill = 'white')
p <- p + scale_fill_manual(values = color_palette)
p <- p + xlab('') + ylab('hME') + NoLegend()
ggsave("F3A_2.png", plot = p, width = 8, height = 4, dpi = 300)

png("SF4B_5.png", width = 8 * 600, height = 6 * 600, res = 300)
HubGeneNetworkPlot(
  seurat_obj,
  mods = "all",
  n_hubs = 5,
  n_other = 10,
  edge_prop = 1
)
dev.off()

modules <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3A/modules.csv", header = TRUE, stringsAsFactors = FALSE)
module_counts <- as.data.frame(table(modules$module))
colnames(module_counts) <- c("Module", "Gene_Count")
color_palette <- c(
  'mEPC 08-M1'  = 'purple',
  'mEPC 08-M2'  = 'blue',
  'mEPC 08-M3'  = 'yellow',
  'mEPC 08-M4'  = 'brown',
  'mEPC 08-M5'  = 'turquoise',
  'mEPC 08-M6'  = 'pink',
  'mEPC 08-M7'  = 'green',
  'mEPC 08-M8'  = 'red',
  'mEPC 08-M9'  = 'magenta',
  'mEPC 08-M10' = 'black',
  'mEPC 08-M11' = 'salmon',
  'mEPC 08-M12' = 'tan',
  'mEPC 08-M13' = 'greenyellow',
  'mEPC 08-M14' = 'cyan'
)

module_counts <- module_counts[module_counts$Module %in% names(color_palette), ]

p <- ggplot(module_counts, aes(x = Gene_Count, y = reorder(Module, Gene_Count), fill = Module)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  labs(x = "Gene Count", y = "Module", title = "Gene Count per Module") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "bold")
  )

ggsave("SF4B_6.png", plot = p, width = 4, height = 8, dpi = 300)

message("条形图已保存为 SF4B_6.png")

saveRDS(seurat_obj, file = 'MNFHC_mEPC_hdWGCNA.rds')


## WGCNA
# TIDE
setwd("C:/R/C_Project/debatch_data/intersect")

data <- read.csv("debatch_mRNA_intersect.csv", header = TRUE, sep = ",")

expression_data <- as.matrix(data[, -1])

rownames(expression_data) <- data[, 1]

library(limma)

exprSet <- normalizeBetweenArrays(expression_data)

min_val <- apply(exprSet, 2, min, na.rm = TRUE)
max_val <- apply(exprSet, 2, max, na.rm = TRUE)
normalized_data <- 2 * ((expression_data - min_val) / (max_val - min_val)) - 1

data[, -1] <- normalized_data

output_path_csv_normalized <- "C:/R/C_Project/SF2/limma_minmax_CDMI.csv"
write.csv(data, file = output_path_csv_normalized, row.names = FALSE)
print(paste("Min-Max归一化后的CSV文件已保存至:", output_path_csv_normalized))

output_path_tab_normalized <- "C:/R/C_Project/SF2/limma_minmax_CDMI.txt"
write.table(data, file = output_path_tab_normalized, sep = "\t", row.names = FALSE, col.names = TRUE)
print(paste("Min-Max归一化后的TAB分隔符文件已保存至:", output_path_tab_normalized))

library(data.table)

file_path <- "C:/R/C_Project/SF2/limma_minmax_TIDE_CDMI.csv"
cdmi_file_path <- "C:/R/C_Project/SF2/CDMI.csv"
output_file_path <- "C:/R/C_Project/SF2/limma_minmax_TIDE_CDMI_reordered.csv"

df <- fread(file_path)

df <- df[complete.cases(df), ]

df[df == "FALSE"] <- 0
df[df == "TRUE"] <- 1

sample_names <- df[[1]]
project_names <- names(df)

df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], as.numeric)

cols_to_remove <- sapply(df[, 2:ncol(df), with = FALSE], function(x) length(unique(x)) == 1)

df <- df[, c(TRUE, !cols_to_remove), with = FALSE]

df[[1]] <- sample_names

project_names <- project_names[c(TRUE, !cols_to_remove)]

df_transposed <- t(df)

df_transposed <- as.data.frame(df_transposed)
colnames(df_transposed) <- df_transposed[1, ]
df_transposed <- df_transposed[-1, ]
df_transposed <- cbind(Calculation_Project = project_names[-1], df_transposed)

df_transposed <- as.data.table(df_transposed)

setnames(df_transposed, old = colnames(df_transposed)[1], new = "")

cdmi_df <- fread(cdmi_file_path)
sample_order <- colnames(cdmi_df)[-1]

reordered_final_df <- df_transposed[, c(1, match(sample_order, colnames(df_transposed)[-1]) + 1), with = FALSE]

fwrite(reordered_final_df, output_file_path, row.names = FALSE)

print(reordered_final_df)

#EMT-ssGSVA
library(clusterProfiler)
library(GSVA)
library(data.table)

file_path <- "C:/R/C_Project/SF2/EMT_SCORE_FILES/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2023.2.Hs.gmt"
lines <- readLines(file_path)
if (length(lines) > 0 && !nzchar(tail(lines, 1))) {
  writeLines(lines, file_path)
} else {
  writeLines(c(lines, ""), file_path)
}

genesets <- clusterProfiler::read.gmt(file_path)
genesets4gsva <- split(genesets$gene, genesets$term)

expr_file <- "C:/R/C_Project/SF2/EMT_SCORE_FILES/limma_log2_CDMI.csv"
expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
rownames(expr) <- expr[, 1]
expr <- as.matrix(expr[,-1])

gsvaP <- ssgseaParam(
  exprData = expr,
  geneSets = genesets4gsva,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)

gsva_data <- gsva(gsvaP)
expr_geneset <- as.data.frame(gsva_data)

output_file <- "C:/R/C_Project/SF2/EMT_SCORE_FILES/EMT_SCORE_limma_log2_CDMI.csv"
write.csv(expr_geneset, output_file, row.names = TRUE)

print(dim(expr_geneset))
print(head(expr_geneset[, 1:4]))

# WGCNA
if (!require(WGCNA)) {
    BiocManager::install("WGCNA")
}
library(WGCNA)
options(stringsAsFactors = FALSE)
femData <- read.csv("C:/R/C_Project/SF2/CDMI.csv")
gene_names <- femData[, 1]
sample_names <- colnames(femData)[-1]
datExpr0 <- as.data.frame(t(sapply(femData[, -1], as.numeric)))
colnames(datExpr0) <- gene_names
rownames(datExpr0) <- sample_names
anyNA(datExpr0)
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples) > 0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}
output_path <- "C:/R/C_Project/SF2/CDMI_For_WGCNA.csv"
write.csv(datExpr0, file = output_path, row.names = TRUE)
sampleTree <- hclust(dist(datExpr0), method = "average")
plot_output_path <- "C:/R/C_Project/SF2/sampleClustering.png"
png(file = plot_output_path, width = 1200, height = 250)
par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels = FALSE)
abline(h = 100, col = "red")
dev.off()
clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
traitData <- read.csv("C:/R/C_Project/SF2/EMT_PHENOTYPES.csv")
dim(traitData)
names(traitData)
allTraits <- traitData[, c(2:11, 19)]
dim(allTraits)
names(allTraits)
femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, traitData$X)
datTraits <- allTraits[traitRows, ]
rownames(datTraits) <- traitData$X[traitRows]
anyDuplicated(rownames(datTraits))
collectGarbage()
sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
output_file <- "C:/R/C_Project/SF2/sample_dendrogram_and_trait_heatmap.png"
png(filename = output_file, width = 1200, height = 400)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap", dendroLabels = FALSE)
dev.off()
cat("Plot saved to", output_file, "\n")
powers = c(1:10, seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

tiff(filename="C:/R/C_Project/SF2/scale_independence.tiff", width=900, height=500)
par(mfrow = c(1, 2))
cex1 = 1.5
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (Power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = paste("Scale Independence"),
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5,
     cex.sub = 1.5)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = cex1,
     col = "red")
abline(h = 0.9, col = "red")
dev.off()

tiff(filename="C:/R/C_Project/SF2/mean_connectivity.tiff", width=500, height=500)
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (Power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean Connectivity"),
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5,
     cex.sub = 1.5)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1,
     col = "red")
dev.off()

softPower = 6

adjacency = adjacency(datExpr, power = softPower)

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

output_file <- "C:/R/C_Project/SF2/gene_clustering_dendrogram.tiff"
tiff(filename = output_file, width = 1200, height = 900)
plot(geneTree,
     xlab = "",
     sub = "",
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)
dev.off()
cat("Plot saved to", output_file, "\n")

minModuleSize = 30

dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

output_file <- "C:/R/C_Project/SF2/gene_dendrogram_and_module_colors.tiff"
tiff(filename = output_file, width = 800, height = 600)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
cat("Plot saved to", output_file, "\n")

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1 - cor(MEs)

METree = hclust(as.dist(MEDiss), method = "average")

output_file <- "C:/R/C_Project/SF2/module_eigengene_clustering.tiff"
tiff(filename = output_file, width = 700, height = 600)
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
MEDissThres = 0.25
abline(h = MEDissThres, col = "red")
dev.off()
cat("Plot saved to", output_file, "\n")

merge = mergeCloseModules(datExpr,
                          dynamicColors,
                          cutHeight = MEDissThres,
                          verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

output_file <- "C:/R/C_Project/SF2/gene_dendrogram_and_module_colors_merged.tiff"
tiff(filename = output_file, width = 1200, height = 740)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Cluster Dendrogram")
dev.off()
cat("Plot saved to", output_file, "\n")
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
output_file <- "C:/R/C_Project/SF2/module_trait_relationships.tiff"
tiff(filename = output_file, width = 1200, height = 900)
par(mar = c(6, 10, 3, 3), cex.lab = 2, cex.axis = 2, cex.main = 2)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1, 1),
               main = paste("Module-Trait Relationships"))
dev.off()
cat("Plot saved to", output_file, "\n")
print(names(datTraits))
if ("EMT_SCORE" %in% names(datTraits)) {
    EMT_SCORE = as.data.frame(datTraits$EMT_SCORE)
    names(EMT_SCORE) = "EMT_SCORE"
} else {
    stop("列 'EMT_SCORE' 在 datTraits 中不存在，请检查列名是否正确。")
}
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
geneTraitSignificance = as.data.frame(cor(datExpr, EMT_SCORE, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(EMT_SCORE), sep = "")
names(GSPvalue) = paste("p.GS.", names(EMT_SCORE), sep = "")
modules = c("black", "yellow")
colors = c("black", "yellow")
moduleGenes = moduleColors %in% modules
output_file <- "C:/R/C_Project/SF2/module_membership_vs_gene_significance_black_yellow.tiff"
tiff(filename = output_file, width = 1200, height = 900)
par(mfrow = c(1, 1))
blackColumn = match("black", modNames)
blackGenes = moduleColors == "black"
verboseScatterplot(abs(geneModuleMembership[blackGenes, blackColumn]),
                   abs(geneTraitSignificance[blackGenes, 1]),
                   xlab = "Module Membership in black module",
                   ylab = "Gene significance for EMT_SCORE",
                   main = "Module membership vs. gene significance\n",
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                   col = "black", pch = 16)
yellowColumn = match("yellow", modNames)
yellowGenes = moduleColors == "yellow"
points(abs(geneModuleMembership[yellowGenes, yellowColumn]),
       abs(geneTraitSignificance[yellowGenes, 1]),
       col = "yellow", pch = 16)
legend("topright", legend = modules, col = colors, pch = 16)
dev.off()
cat("Plot saved to", output_file, "\n")
probes = names(datExpr)
geneInfo0 = data.frame(
  substanceBXH = probes,
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)
modOrder = order(-abs(cor(MEs, EMT_SCORE, use = "p")))
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(
    geneInfo0,
    geneModuleMembership[, modOrder[mod]],
    MMPvalue[, modOrder[mod]]
  )
  names(geneInfo0) = c(
    oldNames,
    paste("MM.",   modNames[modOrder[mod]], sep = ""),
    paste("p.MM.", modNames[modOrder[mod]], sep = "")
  )
}
names(geneTraitSignificance) = paste("GS.",  names(EMT_SCORE), sep = "")
names(GSPvalue)               = paste("p.GS.", names(EMT_SCORE), sep = "")
geneOrder = order(
  geneInfo0$moduleColor,
  -abs(geneInfo0[, paste("GS.", names(EMT_SCORE), sep = "")])
)
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "C:/R/C_Project/SF2/significant_geneset.csv")

nGenes   = ncol(datExpr)
nSamples = nrow(datExpr)
dissTOM  = 1 - TOMsimilarityFromExpr(datExpr, power = 6)
plotTOM  = dissTOM^7
diag(plotTOM) = NA
tiff(
  filename = "C:/R/C_Project/SF2/network_heatmap_all_genes.tiff",
  width    = 900,
  height   = 900
)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network Heatmap Plot，All Genes")
dev.off()

nSelect = 400
set.seed(10)
select     = sample(nGenes, size = nSelect)
selectTOM  = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
plotDiss  = selectTOM^7
diag(plotDiss) = NA
tiff(
  filename = "C:/R/C_Project/SF2/network_heatmap_selected_genes.tiff",
  width    = 900,
  height   = 900
)
TOMplot(plotDiss, selectTree, selectColors, main = "Network Heatmap Plot，Selected Genes")
dev.off()

MEs       = moduleEigengenes(datExpr, moduleColors)$eigengenes
EMT_SCORE = as.data.frame(datTraits$EMT_SCORE)
names(EMT_SCORE) = "EMT_SCORE"
MET       = orderMEs(cbind(MEs, EMT_SCORE))
output_file_dendrogram <- "C:/R/C_Project/SF2/Eigengene_dendrogram.tiff"
output_file_heatmap    <- "C:/R/C_Project/SF2/Eigengene_adjacency_heatmap.tiff"
tiff(filename = output_file_dendrogram, width = 500, height = 750)
par(cex = 1)
plotEigengeneNetworks(
  MET,
  "Eigengene_dendrogram",
  marDendro    = c(0, 4, 2, 0),
  plotHeatmaps = FALSE
)
dev.off()
tiff(filename = output_file_heatmap, width = 500, height = 750)
par(cex = 1)
plotEigengeneNetworks(
  MET,
  "Eigengene_adjacency_heatmap",
  marHeatmap      = c(3, 4, 2, 2),
  plotDendrograms = FALSE,
  xLabelsAngle    = 90
)
dev.off()
cat("Plots saved to", output_file_dendrogram, "and", output_file_heatmap, "\n")

TOM = TOMsimilarityFromExpr(datExpr, power = 6)
dimnames(TOM) = list(names(datExpr), names(datExpr))
modules     = c("yellow", "black")
probes      = names(datExpr)
inModule    = is.finite(match(moduleColors, modules))
modProbes   = probes[inModule]
modTOM      = TOM[modProbes, modProbes]
dimnames(modTOM) = list(modProbes, modProbes)
edgeFilePath = paste(
  "C:/R/C_Project/SF2/CytoscapeInput-edges-",
  paste(modules, collapse = "-"),
  ".txt",
  sep = ""
)
nodeFilePath = paste(
  "C:/R/C_Project/SF2/CytoscapeInput-nodes-",
  paste(modules, collapse = "-"),
  ".txt",
  sep = ""
)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile     = edgeFilePath,
  nodeFile     = nodeFilePath,
  weighted      = TRUE,
  threshold     = 0.02,
  nodeNames     = modProbes,
  altNodeNames  = modProbes,
  nodeAttr      = moduleColors[inModule]
)

##  DEG of paired leision
##  data file is acquired via GEO online tool
file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/Total_Gene_Table.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

data$change <- "N.S."
data$change[data$adj.P.Val < 0.05 & data$logFC >= 0.5] <- "Up"
data$change[data$adj.P.Val < 0.05 & data$logFC <= -0.5] <- "Down"

write.csv(data, file_path, row.names = FALSE)

file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/Total_Gene_Table.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

if (!"Gene.Symbol" %in% colnames(data)) {
  stop("Error: 'Gene.Symbol' 列不存在，请检查数据是否正确。")
}

filtered_data <- data[!grepl("^(RP|MT)", data$Gene.Symbol, ignore.case = TRUE), ]
write.csv(filtered_data, file_path, row.names = FALSE)
cat("已成功删除 Gene.Symbol 以 'RP' 或 'MT' 开头的行，并覆写原文件。\n")

file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/Total_Gene_Table.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

filtered_data <- data[data$change != "N.S.", ]
top_genes <- head(filtered_data[order(filtered_data$logFC, decreasing = FALSE), ], 20)
dat_rep <- data[rownames(data) %in% rownames(top_genes), ]
print(dat_rep)

file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/Total_Gene_Table.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

sig_diff <- data[data$change != "N.S.", ]
cat("筛选后的行数:", nrow(sig_diff), "\n")

volcano_plot <- ggplot(data = data, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  scale_color_manual(values = c("Up" = "seagreen", "N.S." = "darkgray", "Down" = "firebrick3")) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "black", lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.4) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(face = "bold", color = "black", family = "Times", size = 18),
    axis.title.y = element_text(face = "bold", color = "black", family = "Times", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", family = "Times", size = 18),
    plot.subtitle = element_text(hjust = 0.5, family = "Times", size = 12, face = "italic", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black", family = "Times", size = 12)
  ) +
  geom_label_repel(
    data = dat_rep,
    aes(label = Gene.Symbol),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    color = "black",
    show.legend = FALSE
  ) +
  labs(
    x = "log2(Fold_Change)",
    y = "-log10(adj.P.Val)",
    title = "GSE2280 & GSE201777 (Primary_VS_Metastatic)"
  )

output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3C.png"
ggsave(output_path, plot = volcano_plot, width = 8, height = 6.4, dpi = 300)


##  The overlapping part of the candidate gene set
library(UpSetR)

file_paths <- c(
  "hdWGCNA" = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/hdWGCNA.csv",
  "WGCNA"   = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/WGCNA.csv",
  "DEG"     = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/DEG.csv"
)

for (name in names(file_paths)) {
  file_path <- file_paths[[name]]
  if (file.exists(file_path)) {
    cat("\n✅ 读取文件:", name, "\n")
    data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    col_name <- colnames(data)[1]
    value_counts <- table(data[[col_name]])
    num_duplicates <- sum(value_counts > 1)
    cat("列:", col_name, "的重复值个数:", num_duplicates, "\n")
    data_unique <- unique(data)
    write.csv(data_unique, file_path, row.names = FALSE)
    cat("✅ 已删除重复内容并覆写文件:", file_path, "\n")
  } else {
    cat("\n❌ 错误: 文件不存在 ->", file_path, "\n")
  }
}

genes1 <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/hdWGCNA.csv", header = TRUE, stringsAsFactors = FALSE)
genes2 <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/WGCNA.csv", header = TRUE, stringsAsFactors = FALSE)
genes3 <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/DEG.csv", header = TRUE, stringsAsFactors = FALSE)

gene_list1 <- genes1$hdWGCNA
gene_list2 <- genes2$WGCNA
gene_list3 <- genes3$DEG

upsetR_data <- list(
  hdWGCNA = gene_list1,
  WGCNA = gene_list2,
  DEG = gene_list3
)

set_colors <- c("turquoise", "darkgoldenrod1", "firebrick3")

png(filename = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3D.png", width = 1000, height = 800, res = 150)
upset(
  fromList(upsetR_data),
  order.by = "freq",
  main.bar.color = "gray20",
  matrix.color = "skyblue",
  sets.bar.color = set_colors,
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2),
  set_size.angles = -30,
  mainbar.y.label = "Gene Intersections Size"
)
dev.off()

file_paths <- list(
  hdWGCNA = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/hdWGCNA.csv",
  WGCNA   = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/WGCNA.csv",
  DEG     = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/DEG.csv"
)

genes1 <- read.csv(file_paths$hdWGCNA, header = TRUE, stringsAsFactors = FALSE)$hdWGCNA
genes2 <- read.csv(file_paths$WGCNA, header = TRUE, stringsAsFactors = FALSE)$WGCNA
genes3 <- read.csv(file_paths$DEG, header = TRUE, stringsAsFactors = FALSE)$DEG

common_genes <- Reduce(intersect, list(genes1, genes2, genes3))

output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/common_genes.csv"

if (length(common_genes) > 0) {
  write.csv(data.frame(Gene = common_genes), output_file, row.names = FALSE)
  cat("\n✅ 共同基因已保存至:", output_file, "\n")
} else {
  cat("\n⚠️ 没有找到共同基因，未生成 CSV 文件。\n")
}

##  ML part
##  lasso-cox
cdmi_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/CDMI_LASSO_COX.csv"
common_genes_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3C/common_genes.csv"

if (file.exists(common_genes_file)) {
  common_genes <- read.csv(common_genes_file, header = TRUE, stringsAsFactors = FALSE)[[1]]
} else {
  stop("❌ 错误: common_genes.csv 文件不存在！")
}

if (file.exists(cdmi_file)) {
  cdmi_data <- read.csv(cdmi_file, header = TRUE, stringsAsFactors = FALSE)
} else {
  stop("❌ 错误: CDMI_LASSO_COX.csv 文件不存在！")
}

selected_columns <- colnames(cdmi_data) %in% c("X", "OS", "OS_time", common_genes)
filtered_data <- cdmi_data[, selected_columns, drop = FALSE]
write.csv(filtered_data, cdmi_file, row.names = FALSE)

cat("\n✅ 处理完成，已覆盖源文件:", cdmi_file, "\n")

library(dplyr)
library(glmnet)
library(caret)
library(ggplot2)
library(survival)
library(survminer)
library(Matrix)
library(sampling)

file_path <- cdmi_file
data <- read.csv(file_path)

head(data)
table(data$OS)

na_rows <- data[!complete.cases(data$OS_time), ]
if (nrow(na_rows) > 0) {
  print(na_rows)
}

data <- na.omit(data)
head(data)
write.csv(data, file_path, row.names = FALSE)

set.seed(123)
train_id <- strata(data, "OS", size = round(table(data$OS) * 0.62))$ID_unit
trainData <- data[train_id, ]
testData <- data[-train_id, ]
head(trainData)
head(testData)

train_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/CDMI_LASSO_COX_train.csv"
test_path  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/CDMI_LASSO_COX_test.csv"
write.csv(trainData, train_path, row.names = FALSE)
write.csv(testData, test_path, row.names = FALSE)

train.x <- as.matrix(trainData[, 2:33])
train.y <- Surv(trainData$OS_time, trainData$OS)
head(train.x)
head(train.y)

set.seed(123)
train.x <- scale(train.x)
fit <- glmnet(train.x, train.y, family = "cox", alpha = 1)

plot_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_1.png"
png(plot_path, width = 2400, height = 1600, res = 300)
plot(fit, label = FALSE, xvar = "lambda")
dev.off()

set.seed(123)
train.x <- scale(train.x)
fit_cv <- cv.glmnet(train.x, train.y, family = "cox", alpha = 1)

plot_path_cv <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_2.png"
png(plot_path_cv, width = 2400, height = 1600, res = 300)
plot(fit_cv)
dev.off()

print(fit_cv$lambda.min)
print(fit_cv$glmnet.fit)

coef_matrix <- as.matrix(coef(fit_cv, s = "lambda.min"))
print(coef_matrix)

non_zero_features    <- coef_matrix[coef_matrix != 0, , drop = FALSE]
feature_names        <- rownames(non_zero_features)
feature_coefficients <- as.numeric(non_zero_features)
feature_df           <- data.frame(Feature = feature_names, Coefficient = feature_coefficients)
write.csv(feature_df, "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/GENE_LASSO_COX.csv", row.names = FALSE)

glmnet_fit <- glmnet(train.x, train.y, family = "cox", alpha = 1, lambda = fit_cv$lambda.min)
print(glmnet_fit)

coxph_fit <- coxph(train.y ~ train.x)

plot_path_f3 <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_3.png"
png(plot_path_f3, width = 1440, height = 1920, res = 300)
par(cex = 0.75)
plot(
  coef(glmnet_fit),
  coef(coxph_fit),
  col = "red",
  cex = 2,
  xlab = "LASSO Cox Coefficients",
  ylab = "Cox PH Coefficients",
  main = "Coefficient Comparison"
)
abline(0, 1, lty = 2, lwd = 2)
dev.off()

f <- survfit(glmnet_fit, s = 0.05, x = train.x, y = train.y)
survival_plot_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_4.png"
png(survival_plot_path, width = 1600, height = 960, res = 300)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(
  f,
  col = "red",
  lwd = 2,
  ylab = "Survival Probability",
  xlab = "Time",
  cex.axis = 0.75,
  cex.lab = 0.75,
  cex.main = 0.75
)
dev.off()
library(glmnet)
library(survival)

train.surv <- Surv(trainData$OS_time, trainData$OS)
pred <- predict(glmnet_fit, newx = train.x, type = "link")
lasso_cindex <- concordance(train.surv ~ pred)$concordance
print(lasso_cindex)
cox_cindex <- summary(coxph_fit)$concordance[1]
print(cox_cindex)
cindex_results <- data.frame(
  Model = c("LASSO Cox", "Cox PH"),
  Cindex = c(lasso_cindex, cox_cindex)
)
write.csv(
  cindex_results,
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/cindex_results.csv",
  row.names = FALSE
)

library(glmnet)
library(survival)

test.x <- as.matrix(testData[, 2:33])
test.y <- Surv(testData$OS_time, testData$OS)
results <- assess.glmnet(glmnet_fit, newx = test.x, newy = test.y)
print(results)
write.csv(
  data.frame(Measure = names(results), Value = unlist(results)),
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/assessment_results.csv",
  row.names = FALSE
)

## RSF
library(randomForestSRC)
library(ggRandomForests)
library(ggplot2)

file_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/CDMI_LASSO_COX.csv"
rt <- read.csv(file_path, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

colnames(rt)[colnames(rt) == "OS"] <- "fustat"
colnames(rt)[colnames(rt) == "OS_time"] <- "futime"

rfsrc_pbcmy <- rfsrc(
  Surv(futime, fustat) ~ .,
  data = rt,
  nsplit = 10,
  na.action = "na.impute",
  tree.err = TRUE,
  splitrule = "logrank",
  proximity = TRUE,
  forest = TRUE,
  ntree = 1000,
  importance = TRUE
)
print(rfsrc_pbcmy)

gg_dta <- gg_vimp(rfsrc_pbcmy)

vimp_data <- data.frame(
  Feature = gg_dta$vars,
  Importance = gg_dta$vimp
)
write.csv(vimp_data, "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/GENE_VIMP.csv", row.names = FALSE)

plot_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3E_1.png"
plot_gg <- plot(gg_dta)
ggsave(plot_path, plot_gg, width = 8, height = 6)

gg_dta <- gg_minimal_depth(rfsrc_pbcmy)

plot_gg <- plot(gg_dta) +
  theme(
    axis.text.x = element_text(size = rel(0.5)),
    axis.text.y = element_text(size = rel(0.5))
  )
plot_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3E_2.png"
ggsave(plot_path, plot = plot_gg, width = 8, height = 6)

gg_dta_vimp_depth <- gg_minimal_vimp(rfsrc_pbcmy)
gg_dta_vimp_depth <- gg_dta_vimp_depth[rev(seq_len(nrow(gg_dta_vimp_depth))), ]

plot_gg <- plot(gg_dta_vimp_depth) +
  theme(
    axis.text.x = element_text(size = rel(0.5)),
    axis.text.y = element_text(size = rel(0.5))
  )
plot_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3E_3.png"
ggsave(plot_path, plot = plot_gg, width = 4, height = 6)

output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3E-G/RSF.csv"
rsf_data <- gg_dta_vimp_depth[, c("names", "vimp", "depth")]
colnames(rsf_data) <- c("Gene", "VIMP RANK", "MinimalDepth")
write.csv(rsf_data, output_file, row.names = FALSE)
cat("\n✅ 变量重要性数据已保存至:", output_file, "\n")

files_to_delete <- c(
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_4.png",
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3F_3.png",
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3E_1.png",
  "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/figures/F3E_2.png"
)
for (file in files_to_delete) {
  if (file.exists(file)) {
    file.remove(file)
    cat("✅ 已删除:", file, "\n")
  } else {
    cat("⚠️ 文件不存在，跳过:", file, "\n")
  }
}

## WDR54 UMAP projection
import os
import scanpy as sc
import pandas as pd
import harmonypy as hm
import numpy as np
import matplotlib.pyplot as plt
import shutil
import glob
from matplotlib.colors import Normalize
import scipy.sparse as sp

work_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3H"
if not os.path.exists(work_path):
    os.makedirs(work_path)
    print(f"路径 '{work_path}' 已创建。")
else:
    print(f"路径 '{work_path}' 已存在。")
os.chdir(work_path)
print(f"当前工作目录已切换至: {os.getcwd()}")

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F2/F2C/MNFHC_for_RNAvelo.h5ad"
adata = sc.read(file_path)

filtered_cells = adata.obs[adata.obs['leiden_mEPC_CAF'].str.contains('mEPC')]
filtered_indices = filtered_cells.index

if 'WDR54' not in adata.raw.var_names:
    raise ValueError("基因 'WDR54' 不存在于 adata.raw 中。")

wdr54_expression = adata.raw[filtered_indices, 'WDR54'].X.toarray().flatten()

output_df = pd.DataFrame({
    'cell_index': filtered_indices,
    'WDR54_expression': wdr54_expression
})
output_csv_path = "filtered_cells_WDR54_expression.csv"
output_df.to_csv(output_csv_path, index=False)
print(f"结果已保存到: {output_csv_path}")

h5ad_file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F1/F1BC/MNFHC_MalignantEPC.h5ad"
adata = sc.read(h5ad_file_path)

adata.obs['temp_index'] = adata.obs.index.str.replace(r'^(Dis_|LN_)', '', regex=True)
print(adata.obs[['temp_index']].head())

file_path = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F3/F3H/filtered_cells_WDR54_expression.csv"
csv_data = pd.read_csv(file_path)
csv_data.set_index('cell_index', inplace=True)

adata.obs['WDR54_expression'] = adata.obs['temp_index'].map(csv_data['WDR54_expression'])
print(adata.obs[['temp_index', 'WDR54_expression']].head())

wdr54_index = np.where(adata.var_names == "WDR54")[0][0]
adata.obs['WDR54_expression_normalized'] = np.nan
wdr54_expression_values = adata.obs['WDR54_expression']
if sp.issparse(adata.X):
    adata.obs.loc[wdr54_expression_values != 0, 'WDR54_expression_normalized'] = (
        adata.X[wdr54_expression_values != 0, wdr54_index].toarray().flatten()
    )
else:
    adata.obs.loc[wdr54_expression_values != 0, 'WDR54_expression_normalized'] = (
        adata.X[wdr54_expression_values != 0, wdr54_index]
    )

n_neighbors = 100
resolution = 0.5
n_pcs = 30
leiden_mEPC = f'leidenEPC_n{n_neighbors}_r{resolution:.1f}'

umap_coords = adata.obsm['X_umap']
adata.obs[leiden_mEPC] = adata.obs[leiden_mEPC].astype(str)

df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
df[leiden_mEPC] = adata.obs[leiden_mEPC].values

centers = df.groupby(leiden_mEPC)[['UMAP1', 'UMAP2']].mean().reset_index()
centers = centers.rename(columns={'UMAP1': 'Center_UMAP1', 'UMAP2': 'Center_UMAP2'})

df = df.merge(centers, on=leiden_mEPC, how='left')
df['distance_to_center'] = np.sqrt((df['UMAP1'] - df['Center_UMAP1'])**2 + (df['UMAP2'] - df['Center_UMAP2'])**2)

thresholds = df.groupby(leiden_mEPC)['distance_to_center'].quantile(0.95).reset_index()
thresholds = thresholds.rename(columns={'distance_to_center': 'threshold'})

df = df.merge(thresholds[[leiden_mEPC, 'threshold']], on=leiden_mEPC, how='left')
df['is_outlier'] = df['distance_to_center'] > df['threshold']
adata.obs['is_outlier'] = df['is_outlier'].values

outlier_count = adata.obs['is_outlier'].sum()
total_count = adata.n_obs
print(f"检测到 {outlier_count} 个离群点，占总数的 {outlier_count / total_count * 100:.2f}%")

adata_filtered = adata[~adata.obs['is_outlier']].copy()

valid_values = adata_filtered.obs['WDR54_expression_normalized'].dropna()
lower_bound = np.percentile(valid_values, 1)
upper_bound = np.percentile(valid_values, 99)

cmap = plt.cm.viridis

sc.pl.umap(
    adata_filtered,
    color='WDR54_expression_normalized',
    title='UMAP Colored by WDR54 Expression',
    save='F3H_1.png',
    show=False,
    size=30,
    alpha=0.4,
    legend_loc=None,
    cmap=cmap,
    vmin=lower_bound,
    vmax=upper_bound
)

sc.pl.umap(
    adata_filtered,
    color='WDR54_expression_normalized',
    title='UMAP Colored by WDR54 Expression',
    save='F3H_2.png',
    show=False,
    size=30,
    alpha=0.4,
    legend_loc=None,
    cmap=cmap,
    vmin=lower_bound,
    vmax=upper_bound,
    colorbar_loc=None
)
## relation between WDR54 and phenotype
file_path <- "C:/R/C_Project/F2/relation/WDR54_Phenotype.csv"
wdr54_data <- read.csv(file_path, row.names = 1)
wdr54_expr <- wdr54_data$WDR54
sample_names <- rownames(wdr54_data)
min_val <- min(wdr54_expr)
max_val <- max(wdr54_expr)
interval <- (max_val - min_val) / 3
low_threshold <- min_val + interval
medium_threshold <- min_val + 2 * interval
colors <- ifelse(
  wdr54_expr <= low_threshold,
  rgb(255, 224, 193, maxColorValue = 255),
  ifelse(
    wdr54_expr <= medium_threshold,
    rgb(254, 160, 64, maxColorValue = 255),
    rgb(255, 97, 0, maxColorValue = 255)
  )
)
output_file <- "C:/R/C_Project/F2/relation/WDR54_expression_barplot_all_samples_colored.png"
png(
  filename = output_file,
  width = 12000 * 3 / 4 * 1.3,
  height = 600 * 3 / 2 * 1.5 * 1.5
)
bar_width <- 0.5 / length(wdr54_expr) * 1.2 * 1.3
font_size <- 4.0
barplot(
  wdr54_expr,
  main = "WDR54 Expression Levels (All Samples)",
  xlab = "Samples",
  ylab = "Expression Level",
  col = colors,
  border = "white",
  las = 2,
  cex.names = font_size,
  cex.lab = font_size,
  cex.axis = font_size,
  cex.main = font_size,
  names.arg = rep("", length(wdr54_expr)),
  space = -0.3,
  width = bar_width
)
axis(2, lwd = 4, cex.axis = font_size)
legend(
  "topright",
  legend = c("Low", "Medium", "High"),
  fill = c(
    rgb(255, 224, 193, maxColorValue = 255),
    rgb(254, 160, 64,  maxColorValue = 255),
    rgb(255, 97,  0,   maxColorValue = 255)
  ),
  cex = font_size * 0.4,
  bty = "n"
)
dev.off()

library(pheatmap)
file_path_wdr54 <- "C:/R/C_Project/F2/relation/WDR54_Phenotype.csv"
wdr54_data <- read.csv(file_path_wdr54, row.names = 1)
heatmap_data_wdr54 <- t(wdr54_data[, 1:5])
wdr54_data_sorted <- wdr54_data[order(wdr54_data$WDR54), ]
heatmap_data_wdr54_sorted <- t(wdr54_data_sorted[, 1:5])
output_file_wdr54 <- "C:/R/C_Project/F2/relation/WDR54_Heatmap.png"
png(filename = output_file_wdr54, width = 1600, height = 800 / 4)
pheatmap(
  heatmap_data_wdr54_sorted,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend = TRUE,
  fontsize = 12,
  annotation_legend_side = "left",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
dev.off()

library(pheatmap)
file_path <- "C:/R/C_Project/F2/relation/SPP1_Phenotype.csv"
SPP1_data <- read.csv(file_path, row.names = 1)
heatmap_data <- t(SPP1_data[, 1:5])
SPP1_data_sorted <- SPP1_data[order(SPP1_data$SPP1), ]
heatmap_data_sorted <- t(SPP1_data_sorted[, 1:5])
output_file <- "C:/R/C_Project/F2/relation/SPP1_Heatmap.png"
png(filename = output_file, width = 1600, height = 800 / 4)
pheatmap(
  heatmap_data_sorted,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend = TRUE,
  fontsize = 12,
  annotation_legend_side = "left",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
dev.off()
