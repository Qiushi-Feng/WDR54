## KEGG & GO for WGCNA
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

diff <- read.csv(file="C:/R/C_Project/F1/genelist_yellow_and_black.csv")
colnames(diff)[colnames(diff) == 'Gene'] <- 'SYMBOL'

gene.df <- bitr(diff$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene.df)

output_path <- "C:/R/C_Project/F1/gene_conversion.csv"
write.csv(gene.df, file = output_path, row.names = FALSE)
file.exists(output_path)

gene <- gene.df$ENTREZID

ego_ALL <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_CC  <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",  pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_BP  <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",  pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_MF  <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",  pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

output_path <- "C:/R/C_Project/F1"
ego_ALL        <- as.data.frame(ego_ALL)
ego_result_BP  <- as.data.frame(ego_BP)
ego_result_CC  <- as.data.frame(ego_CC)
ego_result_MF  <- as.data.frame(ego_MF)
ego            <- rbind(ego_result_BP, ego_result_CC, ego_result_MF)

write.csv(ego_ALL, file = file.path(output_path, "ego_ALL.csv"),       row.names = TRUE)
write.csv(ego_result_BP, file = file.path(output_path, "ego_result_BP.csv"), row.names = TRUE)
write.csv(ego_result_CC, file = file.path(output_path, "ego_result_CC.csv"), row.names = TRUE)
write.csv(ego_result_MF, file = file.path(output_path, "ego_result_MF.csv"), row.names = TRUE)
write.csv(ego, file = file.path(output_path, "ego.csv"),              row.names = TRUE)

display_number <- c(10, 10, 10)
ego_result_BP <- ego_result_BP[1:display_number[1], ]
ego_result_CC <- ego_result_CC[1:display_number[2], ]
ego_result_MF <- ego_result_MF[1:display_number[3], ]

go_enrich_df <- data.frame(
  ID          = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
  GeneNumber  = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type        = factor(
                  c(
                    rep("biological process", display_number[1]),
                    rep("cellular component", display_number[2]),
                    rep("molecular function", display_number[3])
                  ),
                  levels = c("biological process", "cellular component", "molecular function")
                )
)

for (i in seq_len(nrow(go_enrich_df))) {
  split_desc <- strsplit(go_enrich_df$Description[i], split = " ")[[1]]
  go_enrich_df$Description[i] <- paste(split_desc[1:5], collapse = " ")
}
go_enrich_df$Description <- gsub("NA", "", go_enrich_df$Description)

go_enrich_df$type_order <- factor(
  rev(seq_len(nrow(go_enrich_df))),
  labels = rev(go_enrich_df$Description)
)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

p1 <- ggplot(go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values=COLS) +
  coord_flip() +
  xlab("GO term") + ylab("Gene_Number") +
  labs(title="The Most Enriched GO Terms") +
  theme_bw()
ggsave("C:/R/C_Project/F1/GO_horizontal_barplot.png", plot=p1, width=10, height=8)

p2 <- ggplot(go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values=COLS) +
  theme_bw() +
  xlab("GO term") + ylab("Num of Genes") +
  labs(title="The Most Enriched GO Terms") +
  theme(
    axis.text.x = element_text(face="bold", color="gray50", angle=70, vjust=1, hjust=1)
  )
ggsave("C:/R/C_Project/F1/GO_vertical_barplot.png", plot=p2, width=10, height=8)

kk <- enrichKEGG(gene=gene, keyType="kegg", organism="human", qvalueCutoff=0.05, pvalueCutoff=0.05)
hh <- as.data.frame(kk)
rownames(hh) <- seq_len(nrow(hh))
hh$order <- factor(
  rev(seq_len(nrow(hh))),
  labels = rev(hh$Description)
)

p1 <- ggplot(hh, aes(y=order, x=Count, fill=p.adjust)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_gradient(low="red", high="blue") +
  labs(title="KEGG Pathways Enrichment", x="Gene numbers", y="Pathways") +
  theme(
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    legend.title = element_text(face="bold", size=16)
  ) +
  theme_bw()
ggsave("C:/R/C_Project/F1/KEGG_barplot.png", plot=p1, width=10, height=8)

hh <- as.data.frame(kk)
rownames(hh) <- seq_len(nrow(hh))
hh$order <- factor(
  rev(seq_len(nrow(hh))),
  labels = rev(hh$Description)
)

p2 <- ggplot(hh, aes(y=order, x=Count)) +
  geom_point(aes(size=Count, color=-p.adjust)) +
  scale_color_gradient(low="green", high="red") +
  labs(color=expression(p.adjust), size="Count", x="Gene Number", y="Pathways", title="KEGG Pathway Enrichment") +
  theme_bw()
ggsave("C:/R/C_Project/F1/KEGG_bubbleplot.png", plot=p2, width=10, height=8)

if (!requireNamespace("enrichplot", quietly=TRUE)) install.packages("enrichplot")
library(enrichplot)

kk_sim <- pairwise_termsim(kk)
emap   <- emapplot(kk_sim, showCategory=30)
ggsave("C:/R/C_Project/F1/KEGG_emapplot.png", plot=emap, width=10, height=8)

cnet <- cnetplot(kk, showCategory=5)
ggsave("C:/R/C_Project/F1/KEGG_cnetplot.png", plot=cnet, width=10, height=8)


## GSEA(GO & Hallmarker)
#——————————————————————————————
# step1. 差异基因分析
#——————————————————————————————
library(limma)
library(edgeR)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(viridis)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

if (!requireNamespace("enrichplot", quietly = TRUE)) {
    BiocManager::install("enrichplot")
}
library(enrichplot)

if (!requireNamespace("DOSE", quietly = TRUE)) {
    BiocManager::install("DOSE")
}
library(DOSE)

if (!requireNamespace("gseaPlotting", quietly = TRUE)) {
    BiocManager::install("gseaPlotting")
}
library(gseaPlotting)

if ("gseadist" %in% ls("package:gseaPlotting")) {
    cat("gseadist 函数已成功加载并可用。\n")
} else {
    cat("未找到 gseadist 函数。请确认该函数所在的包。\n")
}

set.seed(456)

data <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI.csv", row.names = 1, check.names = FALSE)
WDR54_expression <- as.numeric(data["WDR54", ])
median_value <- median(WDR54_expression)
group <- ifelse(WDR54_expression >= median_value, "High", "Low")
data_with_group <- rbind(data, Group = group)
write.csv(data_with_group, "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_Group.csv", row.names = TRUE)

setwd("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA")
cdmi_data <- read.csv("CDMI.csv", row.names = 1)
cdmi_genes <- rownames(cdmi_data)
gpl_data <- read.csv("GPL10558-50081.csv")
gpl_symbols <- gpl_data$Symbol
gpl_entrez_ids <- gpl_data$Entrez_Gene_ID
matched_indices <- which(gpl_symbols %in% cdmi_genes)
matched_symbols <- gpl_symbols[matched_indices]
matched_entrez_ids <- gpl_entrez_ids[matched_indices]
output_data <- data.frame(Symbol = matched_symbols, Entrez_Gene_ID = matched_entrez_ids)
output_data <- output_data[!duplicated(output_data$Symbol), ]
write.csv(output_data, file = "Entrez_Reflection.csv", row.names = FALSE)

unmatched_genes <- setdiff(cdmi_genes, gpl_symbols)
cat("无匹配的基因名如下：\n")
print(unmatched_genes)

data <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_Group.csv", row.names = 1)

group <- as.factor(as.character(data["Group", ]))
levels(group) <- c("Low", "High")

exprSet <- data[-which(rownames(data) == "Group"), ]
gene_names <- rownames(exprSet)
exprSet <- as.data.frame(lapply(exprSet, as.numeric))
rownames(exprSet) <- gene_names

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cat("设计矩阵前几行：\n")
print(head(design))

fit <- lmFit(exprSet, design)
contrast.matrix <- makeContrasts(High-Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

diff_limma <- topTable(fit2, adjust = "fdr", number = Inf)
cat("差异表达分析结果行名：\n")
print(rownames(diff_limma))

write.csv(diff_limma, file = "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_diff_limma.csv", row.names = TRUE)
head(diff_limma)


#——————————————————————————————
# step2. C5GSEA计算以及通路初筛与山脊图绘制
#——————————————————————————————
diff_limma <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_diff_limma.csv", row.names = 1)

gene_entrezid <- bitr(geneID = rownames(diff_limma),
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

gene_entrezid <- merge(gene_entrezid, diff_limma, by.x = "SYMBOL", by.y = "row.names")

genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene)

gsea_res <- GSEA(genelist,
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 seed = 456)

gsea_res_df <- as.data.frame(gsea_res)

out_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA"
out_file <- file.path(out_dir, "CDMI_GO_GSEA_results.csv")

write.csv(gsea_res_df, file = out_file, row.names = TRUE)
message("已将 GSEA 结果保存至：", out_file)

head(gsea_res_df)


## 将GSEA的结果中的ID输入GPT，挑选出和自己工作相关的路径保存整理到CDMI_GSEA_results_Selected.csv（这个工作量非常庞大）

### 通路挑选策略
#0.策略制定：分为迁移、EMT、免疫三个方面，在正图中放置同向改变的前三重要的项目，附图中放确凿的气泡图
#1. 将所有富集到的ID通过GPT翻译并保存到 GSEA结果.doc
#2. 挑选其中与EMT、迁移、免疫有关的词条，将对应数据分别摘抄到
#CDMI_GSEA_results_EMT.csv/CDMI_GSEA_results_Migration.csv/CDMI_GSEA_results_Immune.csv
#3. 删除其中不符合预测变化趋势的词条，以及较为边缘的，功能模糊的词条
#4. 各挑选重要程度高，变化方向一致的五个词条，放入CDMI_GSEA_results_EMT_Curve.csv/
#CDMI_GSEA_results_Migration_Curve.csv/CDMI_GSEA_results_Immune_Curve.csv文件中

selected_pathways <- c(
  "HP_SQUAMOUS_CELL_CARCINOMA",
  "HP_NEOPLASM_OF_HEAD_AND_NECK",
  "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION",
  "GOBP_MESENCHYMAL_CELL_PROLIFERATION",
  "GOBP_MESENCHYME_DEVELOPMENT",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
  "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_SUPERFAMILY_BINDING",
  "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_BINDING"
)

gsea_res@result$Description <- factor(
  gsea_res@result$Description,
  levels = rev(selected_pathways)
)

ids <- gsea_res@result$ID[match(selected_pathways, gsea_res@result$Description)]

ridgeplot_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_EMT.PNG"
png(filename = ridgeplot_output_file, width = 600, height = 800)
ridgeplot(
  gsea_res,
  showCategory = selected_pathways,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 30,
  orderBy = "Description",
  decreasing = FALSE
) + theme(axis.text.y = element_text(size = 8))
dev.off()

density_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_EMT_density_plot.PNG"
png(filename = density_output_file, width = 600, height = 800)
gseadist(
  gsea_res,
  IDs = ids,
  type = "density"
) + theme(
  legend.direction = "vertical",
  legend.text = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  line = element_line(linewidth = 1.2)
)
dev.off()

selected_pathways <- c(
  "GOBP_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_REGULATION_OF_T_HELPER_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_B_CELL_DIFFERENTIATION",
  "GOBP_B_CELL_DIFFERENTIATION",
  "GOBP_HUMORAL_IMMUNE_RESPONSE",
  "GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY"
)

gsea_res@result$Description <- factor(
  gsea_res@result$Description,
  levels = rev(selected_pathways)
)

ids <- gsea_res@result$ID[match(selected_pathways, gsea_res@result$Description)]
ridgeplot_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Immune.PNG"
png(filename = ridgeplot_output_file, width = 600, height = 800)
ridgeplot(
  gsea_res,
  showCategory = selected_pathways,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 30,
  orderBy = "Description",
  decreasing = FALSE
) +
  scale_fill_gradient(
    low  = rgb(195,   0, 120, maxColorValue = 255),
    high = rgb( 84, 159, 154, maxColorValue = 255)
  ) +
  theme(axis.text.y = element_text(size = 8))
dev.off()

density_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Immune_density_plot.PNG"
png(filename = density_output_file, width = 600, height = 800)
gseadist(
  gsea_res,
  IDs  = ids,
  type = "density"
) +
  theme(
    legend.direction = "vertical",
    legend.text      = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    line             = element_line(linewidth = 1.2)
  )
dev.off()

selected_pathways <- c(
  "GOBP_CELL_MIGRATION",
  "GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",
  "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY",
  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
  "GOBP_COLLAGEN_METABOLIC_PROCESS",
  "GOBP_COLLAGEN_CATABOLIC_PROCESS",
  "GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS",
  "GOBP_REGULATION_OF_CYTOSKELETON_ORGANIZATION",
  "GOMF_METALLOPEPTIDASE_ACTIVITY"
)

gsea_res@result$Description <- factor(
  gsea_res@result$Description,
  levels = rev(selected_pathways)
)
ids <- gsea_res@result$ID[match(selected_pathways, gsea_res@result$Description)]

ridgeplot_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Migration.PNG"
png(filename = ridgeplot_output_file, width = 600, height = 800)
ridgeplot(
  gsea_res,
  showCategory = selected_pathways,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 30,
  orderBy = "Description",
  decreasing = FALSE
) +
  scale_fill_gradient(
    low  = rgb(243, 215,   7, maxColorValue = 255),
    high = rgb( 89, 183, 143, maxColorValue = 255)
  ) +
  theme(axis.text.y = element_text(size = 8))
dev.off()

density_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Migration_density_plot.PNG"
png(filename = density_output_file, width = 600, height = 800)
gseadist(
  gsea_res,
  IDs  = ids,
  type = "density"
) +
  theme(
    legend.direction = "vertical",
    legend.text      = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    line             = element_line(linewidth = 1.2)
  )
dev.off()

selected_pathways <- c(
  "GOBP_REGULATION_OF_T_CELL_MIGRATION",
  "GOBP_NEGATIVE_REGULATION_OF_CHEMOTAXIS",
  "GOBP_POSITIVE_REGULATION_OF_T_CELL_MIGRATION",
  "GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION",
  "GOBP_POSITIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "GOBP_MYELOID_LEUKOCYTE_ACTIVATION"
)

gsea_res@result$Description <- factor(
  gsea_res@result$Description,
  levels = rev(selected_pathways)
)
ids <- gsea_res@result$ID[match(selected_pathways, gsea_res@result$Description)]

ridgeplot_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Immune_Migration.PNG"
png(filename = ridgeplot_output_file, width = 600, height = 800)
ridgeplot(
  gsea_res,
  showCategory = selected_pathways,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 30,
  orderBy = "Description",
  decreasing = FALSE
) +
  scale_fill_gradient(
    low  = rgb(237,  52,  47, maxColorValue = 255),
    high = rgb( 23, 134,  66, maxColorValue = 255)
  ) +
  theme(axis.text.y = element_text(size = 8))
dev.off()

density_output_file <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest/GSEA_Immune_Migration_density_plot.PNG"
png(filename = density_output_file, width = 600, height = 800)
gseadist(
  gsea_res,
  IDs  = ids,
  type = "density"
) +
  theme(
    legend.direction = "vertical",
    legend.text      = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    line             = element_line(linewidth = 1.2)
  )
dev.off()

#——————————————————————————————
# step3. C5基因集GSEA分析，基因集筛选与可视化
#——————————————————————————————
diff_limma <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_diff_limma.csv", row.names = 1)

gene_entrezid <- bitr(
  geneID = rownames(diff_limma),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

gene_entrezid <- merge(
  gene_entrezid,
  diff_limma,
  by.x = "SYMBOL",
  by.y = "row.names"
)

genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene)

gsea_res <- GSEA(
  genelist,
  TERM2GENE = m_t2g,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  seed = 456
)

output_dir <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/pretest"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

gsea_res_symbol <- setReadable(gsea_res, "org.Hs.eg.db", "ENTREZID")

geneSetIDs <- c(
  "HP_SQUAMOUS_CELL_CARCINOMA", "HP_NEOPLASM_OF_HEAD_AND_NECK",
  "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION",
  "GOBP_MESENCHYMAL_CELL_PROLIFERATION", "GOBP_MESENCHYME_DEVELOPMENT",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
  "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_SUPERFAMILY_BINDING",
  "GOMF_TUMOR_NECROSIS_FACTOR_RECEPTOR_BINDING",
  "GOBP_CELL_MIGRATION", "GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",
  "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY", "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
  "GOBP_COLLAGEN_METABOLIC_PROCESS", "GOBP_COLLAGEN_CATABOLIC_PROCESS",
  "GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_CYTOSKELETON_ORGANIZATION",
  "GOMF_METALLOPEPTIDASE_ACTIVITY", "GOBP_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_REGULATION_OF_T_HELPER_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_B_CELL_DIFFERENTIATION", "GOBP_B_CELL_DIFFERENTIATION",
  "GOBP_HUMORAL_IMMUNE_RESPONSE",
  "GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY"
)

for (geneSetID in geneSetIDs) {
  geneSetIndex <- which(gsea_res_symbol$ID == geneSetID)
  p <- gseaplot2(
    gsea_res_symbol,
    geneSetID = geneSetIndex,
    title = gsea_res_symbol$Description[geneSetIndex]
  )
  p[[1]] <- p[[1]] + theme(title = element_text(color = "red"))
  output_path <- file.path(output_dir, paste0(geneSetID, ".png"))
  png(filename = output_path, width = 1000, height = 800)
  print(p)
  dev.off()
}

geneSetIDs <- c(
  "GOBP_CELL_MIGRATION", "GOBP_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_MESENCHYME_DEVELOPMENT"
)
labels <- c("CELL_MIGRATION", "IMMUNE_EFFECTOR", "MESENCHYME_DEVELOPMENT")
custom_colors <- c("#FBB463", "#8DD1C6", "#BD9AAD")
geneSetIndex <- which(gsea_res$ID %in% geneSetIDs)

p <- gseaplot2(gsea_res, geneSetID = geneSetIndex)
p[[1]] <- p[[1]] +
  scale_color_manual(values = custom_colors, labels = labels) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.8, linetype = 2) +
  theme(legend.position = "top")
p[[2]] <- p[[2]] + scale_color_manual(values = custom_colors)
p[[3]] <- p[[3]] +
  geom_hline(yintercept = 0, color = "steelblue", linewidth = 0.5, linetype = 2)

output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/figure/GSEA_GO_Combined_Plot.png"
png(filename = output_path, width = 1500, height = 1200, res = 300)
print(p)
dev.off()
print(p)


#———————————————————————————————————————————————————
# step4. Hallmarker基因集GSEA分析，基因集筛选与可视化
#———————————————————————————————————————————————————
diff_limma <- read.csv("/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/CDMI_diff_limma.csv", row.names = 1)
gene_entrezid <- bitr(
  geneID   = rownames(diff_limma),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = "org.Hs.eg.db"
)
gene_entrezid <- merge(
  gene_entrezid,
  diff_limma,
  by.x = "SYMBOL",
  by.y = "row.names"
)
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist, decreasing = TRUE)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
gsea_res <- GSEA(
  genelist,
  TERM2GENE     = m_t2g,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "none",
  seed          = 456
)
gsea_res_df <- as.data.frame(gsea_res)
out_dir  <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA"
out_file <- file.path(out_dir, "CDMI_hallmarker_GSEA_results.csv")
write.csv(gsea_res_df, file = out_file, row.names = TRUE)
message("已将 GSEA 结果保存至：", out_file)

output_path <- "/media/desk16/tly6105/Essential_Data/Spatial_and_single_cell/F5/GSEA/figure/GSEA_hallmarker_Combined_Plot.png"
geneSetIDs     <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING")
labels         <- c("EMT", "TGF-beta")
custom_colors  <- c("#C22284", "#3A53A4")
geneSetIndex   <- which(gsea_res$ID %in% geneSetIDs)
p <- gseaplot2(gsea_res, geneSetID = geneSetIndex)
p[[1]] <- p[[1]] +
  scale_color_manual(values = custom_colors, labels = labels) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.8, linetype = 2) +
  theme(legend.position = "top")
p[[2]] <- p[[2]] + scale_color_manual(values = custom_colors)
p[[3]] <- p[[3]] + geom_hline(yintercept = 0, color = "steelblue", linewidth = 0.5, linetype = 2)

png(filename = output_path, width = 1500, height = 1200, res = 300)
print(p)
dev.off()
