## WDR54 EXPRESSION
library(ggplot2)
library(ggsignif)
library(ggthemes)
library(reshape2)
library(dplyr)
library(tidyr)
library(limma)

data1 <- read.csv("C:/R/C_Project/F3/WDR54_Contrast/GSE65858_Subtype.csv")
p1 <- ggplot(data1, aes(x = Tissue_type, y = WDR54, fill = Tissue_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = Tissue_type),
                position = position_jitterdodge(jitter.width = 0.03, dodge.width = 0.75),
                size = 1.5) +
    labs(title = "WDR54 Expression Levels Across Tissue Subtypes",
         x = "Tissue Subtype",
         y = "WDR54 Expression Level") +
    scale_fill_manual(values = c("Atypical" = "lightblue", 
                                 "Classical" = "lightgreen", 
                                 "Mesenchymal" = "lightcoral", 
                                 "Basal" = "lightyellow")) +
    scale_color_manual(values = c("Atypical" = "#6888F5", 
                                  "Classical" = "#78C679", 
                                  "Mesenchymal" = "#D77071", 
                                  "Basal" = "#F7D56B")) +
    scale_y_continuous(limits = c(7.75, 11)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgrey"),
          panel.grid.minor = element_line(color = "lightgrey"),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))
ggsave("C:/R/C_Project/F3/WDR54_Contrast/WDR54_Expression_Subtypes_YAxis_7.75_11.png", plot = p1, width = 9.6, height = 6)

data2 <- read.csv("C:/R/C_Project/F3/WDR54_Contrast/CDMI_Tumor_VS_Normal.csv")
p2 <- ggplot(data2, aes(x = GROUP, y = WDR54, fill = GROUP)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = GROUP),
                position = position_jitterdodge(jitter.width = 0.03, dodge.width = 0.75),
                size = 1.5) +
    labs(title = "WDR54 Expression Levels in Normal vs Tumor Samples",
         x = "Sample Group",
         y = "WDR54 Expression Level") +
    scale_fill_manual(values = c("Normal" = "lightblue", "Tumor" = "lightcoral")) +
    scale_color_manual(values = c("Normal" = "#6888F5", "Tumor" = "#D77071")) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgrey"),
          panel.grid.minor = element_line(color = "lightgrey"),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    geom_signif(comparisons = list(c("Normal", "Tumor")),
                map_signif_level = FALSE,
                test = "t.test",
                textsize = 3,
                vjust = 0.2,
                tip_length = 0.02)
ggsave("C:/R/C_Project/F3/WDR54_Contrast/WDR54_Expression_Boxplot_with_PValue.png", plot = p2, width = 5.333, height = 6)

data3 <- read.csv("C:/R/C_Project/F3/WDR54_Contrast/GSE201777_Primary_VS_Metastatic.csv")
p3 <- ggplot(data3, aes(x = Tissue_type, y = WDR54, fill = Tissue_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = Tissue_type),
                position = position_jitterdodge(jitter.width = 0.03, dodge.width = 0.75),
                size = 1.5) +
    labs(title = "WDR54 Expression Levels in Primary vs Metastatic Samples",
         x = "Tissue Type",
         y = "WDR54 Expression Level") +
    scale_fill_manual(values = c("Primary" = "#ACD2C7", "Metastatic" = "#CA8BA8")) +
    scale_color_manual(values = c("Primary" = "#549F9A", "Metastatic" = "#C30078")) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgrey"),
          panel.grid.minor = element_line(color = "lightgrey"),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    geom_signif(comparisons = list(c("Primary", "Metastatic")),
                map_signif_level = FALSE,
                test = "t.test",
                textsize = 3,
                vjust = 0.2,
                tip_length = 0.02)
ggsave("C:/R/C_Project/F3/WDR54_Contrast/WDR54_Expression_Primary_vs_Metastatic.png", plot = p3, width = 5.333, height = 6)

data4 <- read.csv("C:/R/C_Project/F3/WDR54_Contrast/TCGA_Paired_Normal_VS_Tumor.csv")
data_long <- reshape2::melt(data4, id.vars = "Patient_ID", variable.name = "Tissue_Type", value.name = "WDR54_Expression")
stat_test <- t.test(WDR54_Expression ~ Tissue_Type, data = data_long)
p4 <- ggplot(data_long, aes(x = Tissue_Type, y = WDR54_Expression, fill = Tissue_Type)) +
    geom_line(aes(group = Patient_ID), size = 0.5) +
    geom_point(shape = 21, size = 3, stroke = 0.6, color = 'black') +
    scale_fill_manual(values = c("Normal" = "#6888F5", "Tumor" = "#D77071")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text.y = element_text(size = 14, face = "bold", color = "black"),
          axis.text.x = element_text(size = 14, face = "bold", color = "black"),
          axis.title.y = element_text(size = 15, color = "black", face = "bold"),
          axis.title.x = element_text(size = 15, color = "black", face = "bold"),
          axis.line = element_line(color = "black")) +
    labs(x = 'Tissue Type', y = 'WDR54 Expression Level') +
    geom_signif(comparisons = list(c("Normal", "Tumor")),
                map_signif_level = FALSE,
                annotations = paste("p =", signif(stat_test$p.value, digits = 3)),
                y_position = max(data_long$WDR54_Expression) + 0.5,
                tip_length = 0.02)
print(p4)
ggsave("C:/R/C_Project/F3/WDR54_Contrast/WDR54_Expression_Paired_Normal_vs_Tumor_with_Significance.png", plot = p4, width = 5.6, height = 6)



## WDR54 clinical part
library(readr)
input_file_path <- "C:/R/C_Project/F3/Clinical_landscape/TCGA_Clinical_Landscape.csv"
output_file_path <- "C:/R/C_Project/F3/Clinical_landscape/Processed_TCGA_Clinical_Landscape.csv"
clinical_data <- read.csv(input_file_path, row.names = 1, stringsAsFactors = FALSE)
WDR54_values <- as.numeric(clinical_data[,"WDR54"])
WDR54_median <- median(WDR54_values, na.rm = TRUE)
clinical_data[,"WDR54"] <- ifelse(WDR54_values < WDR54_median, "Low", "High")
Age_values <- as.numeric(clinical_data[,"Age"])
clinical_data[,"Age"] <- cut(Age_values, 
                             breaks = c(-Inf, 20, 40, 50, 60, 70, 80, Inf), 
                             labels = c("<20", "[20, 40)", "[40, 50)", "[50, 60)", "[60, 70)", "[70, 80)", "≥80"), 
                             right = FALSE)

Tobacco_exposure_values <- as.numeric(clinical_data[,"Tobacco_exposure"])
clinical_data[,"Tobacco_exposure"] <- ifelse(is.na(Tobacco_exposure_values), 
                                             NA,
                                             as.character(cut(Tobacco_exposure_values, 
                                                              breaks = c(-Inf, 10, 30, 50, Inf), 
                                                              labels = c("<10", "[10, 30)", "[30, 50)", "≥50"), 
                                                              right = FALSE)))
write.csv(clinical_data, file = output_file_path, row.names = TRUE)
cat("处理后的数据已保存至:", output_file_path, "\n")
1-Low [0, 191, 196]
2-High[248 ,118 ,109]
1-Floor of mouth [247, 166, 172]
2-Tongue [238,193,134]
3-Pharynx and larynx [178,219,185]
4-Gum [124,124,186]
5-Tonsil [76,108,67]
6-Other parts of mouth [204,204,204]
1-<40 [255, 192 128]
2-[40, 60) [242, 128, 128]
3-[60, 80) [107, 126, 185]
4-≥80 [114, 195, 163]
1-Female [67, 163, 239]
2-Male [239, 118, 123]
1-Yes [9, 154, 99]
2-No [173, 7, 227]
3-Not Reported [204,204,204]
1-<10 [191, 145, 213]
2-[10, 30) [128, 163, 213]
3-[30, 50) [255, 128, 128]
4-≥50 [128, 208, 195]
5-NA [204,204,204]
2-Negative [149, 169, 56]
3-Positive [173. 54, 136]
4-NA [204,204,204]
1-T1 [236, 210, 211]
2-T2 [219, 131, 130]
3-T3 [203, 24, 28]
4-T4 [108, 14, 15]
5-Tx [204,204,204]
1-N0 [211, 208, 253]
2-N1 [132, 121, 249]
3-N2 [70, 38, 246]
4-N3 [41, 17, 115]
5-Nx [204,204,204]
1-M0 [202, 193,212]
2-M1 [97, 26,68]
3-Mx [204,204,204]
1-i [218, 233, 228]
2-Ii [184, 219, 179]
3-Iii [122, 182, 86]
4-Iv [21, 100, 52]
5-Not Reported [204,204,204]
library(ComplexHeatmap)
library(circlize) 
file_path <- "C:/R/C_Project/F3/Clinical_landscape/primary_site.csv"
clinical_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
wdr54_data <- as.matrix(clinical_data[,"WDR54", drop = FALSE])
my_color <- colorRamp2(c(1, 2), 
                       c(rgb(0, 191, 196, maxColorValue = 255),  
                         rgb(248, 118, 109, maxColorValue = 255)))  
heatmap_plot <- Heatmap(wdr54_data, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE)) 

output_path <- "C:/R/C_Project/F3/Clinical_landscape/WDR54_heatmap_primary_site.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
library(ComplexHeatmap)
library(circlize) 
file_path <- "C:/R/C_Project/F3/Clinical_landscape/primary_site.csv"
clinical_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
primary_site_data <- as.matrix(clinical_data[,"primary_site", drop = FALSE])
my_color <- colorRamp2(c(1, 2, 3, 4, 5, 6), 
                       c(rgb(247, 166, 172, maxColorValue = 255),  
                         rgb(238, 193, 134, maxColorValue = 255),  
                         rgb(178, 219, 185, maxColorValue = 255),  
                         rgb(124, 124, 186, maxColorValue = 255),  
                         rgb(76, 108, 67, maxColorValue = 255),    
                         rgb(204, 204, 204, maxColorValue = 255))) 

heatmap_plot <- Heatmap(primary_site_data, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE))  
output_path <- "C:/R/C_Project/F3/Clinical_landscape/primary_site_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
library(ComplexHeatmap)
library(circlize) 
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Age.csv"
age_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
age_data_matrix <- as.matrix(age_data[,"Age", drop = FALSE])
my_color <- colorRamp2(c(1, 2, 3, 4), 
                       c(rgb(255, 192, 128, maxColorValue = 255), 
                         rgb(242, 128, 128, maxColorValue = 255),
                         rgb(107, 126, 185, maxColorValue = 255),  
                         rgb(114, 195, 163, maxColorValue = 255))) 
heatmap_plot <- Heatmap(age_data_matrix, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE)) 

output_path <- "C:/R/C_Project/F3/Clinical_landscape/Age_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
library(ComplexHeatmap)
library(circlize) 
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Gender.csv"
gender_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
gender_data_matrix <- as.matrix(gender_data[,"Gender", drop = FALSE])
my_color <- colorRamp2(c(1, 2), 
                       c(rgb(67, 163, 239, maxColorValue = 255),  
                         rgb(239, 118, 123, maxColorValue = 255))) 

heatmap_plot <- Heatmap(gender_data_matrix, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE))  

output_path <- "C:/R/C_Project/F3/Clinical_landscape/Gender_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
library(ComplexHeatmap)
library(circlize) 
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Alcohol_exposure.csv"
alcohol_exposure_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
alcohol_exposure_data_matrix <- as.matrix(alcohol_exposure_data[,"Alcohol_exposure", drop = FALSE])
my_color <- colorRamp2(c(1, 2, 3), 
                       c(rgb(9, 154, 99, maxColorValue = 255),  
                         rgb(173, 7, 227, maxColorValue = 255), 
                         rgb(204, 204, 204, maxColorValue = 255))) 
heatmap_plot <- Heatmap(alcohol_exposure_data_matrix, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE))  
output_path <- "C:/R/C_Project/F3/Clinical_landscape/Alcohol_exposure_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Tobacco_exposure.csv"
tobacco_exposure_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

tobacco_exposure_data_matrix <- as.matrix(tobacco_exposure_data[,"Tobacco_exposure", drop = FALSE])
my_color <- colorRamp2(c(1, 2, 3, 4, 5), 
                       c(rgb(191, 145, 213, maxColorValue = 255),  
                         rgb(128, 163, 213, maxColorValue = 255),  
                         rgb(255, 128, 128, maxColorValue = 255),  
                         rgb(128, 208, 195, maxColorValue = 255),  
                         rgb(204, 204, 204, maxColorValue = 255))) 
heatmap_plot <- Heatmap(tobacco_exposure_data_matrix, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE))

output_path <- "C:/R/C_Project/F3/Clinical_landscape/Tobacco_exposure_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
file_path <- "C:/R/C_Project/F3/Clinical_landscape/HPV.csv"
hpv_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
hpv_data_matrix <- as.matrix(hpv_data[,"HPV", drop = FALSE])
my_color <- colorRamp2(c(2, 3, 4), 
                       c(rgb(149, 169, 56, maxColorValue = 255),   
                         rgb(173, 54, 136, maxColorValue = 255),   
                         rgb(204, 204, 204, maxColorValue = 255))) 
heatmap_plot <- Heatmap(hpv_data_matrix, 
                        col = my_color, 
                        show_row_names = FALSE, 
                        show_column_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        heatmap_legend_param = list(title = NULL, show = FALSE)) 
output_path <- "C:/R/C_Project/F3/Clinical_landscape/HPV_heatmap.png"
png(filename = output_path, width = 80, height = 1000)
draw(heatmap_plot)
dev.off()
library(ggplot2)
library(dplyr)
library(reshape2)
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Pathologic_T.csv"
pathologic_t_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
pathologic_t_long <- melt(pathologic_t_data, id.vars = "WDR54", variable.name = "Variable", value.name = "Pathologic_T")
pathologic_t_long <- pathologic_t_long[pathologic_t_long$Variable == "Pathologic_T",]
pathologic_t_long$Pathologic_T <- factor(pathologic_t_long$Pathologic_T, levels = c("1", "2", "3", "4", "5"))
pathologic_t_percentage <- pathologic_t_long %>%
  group_by(WDR54, Pathologic_T) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()
contingency_table <- table(pathologic_t_long$WDR54, pathologic_t_long$Pathologic_T)
chi_test <- chisq.test(contingency_table)
p_value <- chi_test$p.value
cat("P-value for the difference in Pathologic_T distribution between WDR54 groups:", p_value, "\n")
p_value_str <- sprintf("%.3e", p_value)
output_path <- paste0("C:/R/C_Project/F3/Clinical_landscape/WDR54_Pathologic_T_percentage_p_", p_value_str, ".png")
plot <- ggplot(pathologic_t_percentage, aes(x = WDR54, y = Percentage, fill = Pathologic_T)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#FF9999", "#FFCC99", "#99CCFF", "#66CC99", "#CCCCCC")) + 
  labs(title = "Percentage Stacked Bar Chart of Pathologic T by WDR54 Expression",
       x = "WDR54 Expression",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 20, b = 15, l = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 1.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )
original_width <- 8 
new_width <- original_width * 0.618
ggsave(filename = output_path, plot = plot, width = new_width, height = 6, units = "in")
library(ggplot2)
library(dplyr)
library(reshape2)
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Pathologic_N.csv"
pathologic_n_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
pathologic_n_long <- melt(pathologic_n_data, id.vars = "WDR54", variable.name = "Variable", value.name = "Pathologic_N")
pathologic_n_long <- pathologic_n_long[pathologic_n_long$Variable == "Pathologic_N",]
pathologic_n_long$Pathologic_N <- factor(pathologic_n_long$Pathologic_N, levels = c("1", "2", "3", "4", "5"))
pathologic_n_percentage <- pathologic_n_long %>%
  group_by(WDR54, Pathologic_N) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

contingency_table <- table(pathologic_n_long$WDR54, pathologic_n_long$Pathologic_N)
chi_test <- chisq.test(contingency_table)
p_value <- chi_test$p.value
cat("P-value for the difference in Pathologic_N distribution between WDR54 groups:", p_value, "\n")
p_value_str <- sprintf("%.3e", p_value)
output_path <- paste0("C:/R/C_Project/F3/Clinical_landscape/WDR54_Pathologic_N_percentage_p_", p_value_str, ".png")
plot <- ggplot(pathologic_n_percentage, aes(x = WDR54, y = Percentage, fill = Pathologic_N)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#FF9999", "#FFCC99", "#99CCFF", "#66CC99", "#CCCCCC")) + 
  labs(title = "Percentage Stacked Bar Chart of Pathologic N by WDR54 Expression",
       x = "WDR54 Expression",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) +  
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 20, b = 15, l = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 1.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
original_width <- 8
new_width <- original_width * 0.618
ggsave(filename = output_path, plot = plot, width = new_width, height = 6, units = "in")
library(ggplot2)
library(dplyr)
library(reshape2)
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Pathologic_M.csv"
pathologic_m_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
pathologic_m_long <- melt(pathologic_m_data, id.vars = "WDR54", variable.name = "Variable", value.name = "Pathologic_M")
pathologic_m_long <- pathologic_m_long[pathologic_m_long$Variable == "Pathologic_M",]
pathologic_m_long$Pathologic_M <- factor(pathologic_m_long$Pathologic_M, levels = c("1", "2"))
pathologic_m_percentage <- pathologic_m_long %>%
  group_by(WDR54, Pathologic_M) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()
contingency_table <- table(pathologic_m_long$WDR54, pathologic_m_long$Pathologic_M)
chi_test <- chisq.test(contingency_table)
p_value <- chi_test$p.value
cat("P-value for the difference in Pathologic_M distribution between WDR54 groups:", p_value, "\n")
p_value_str <- sprintf("%.3e", p_value)
output_path <- paste0("C:/R/C_Project/F3/Clinical_landscape/WDR54_Pathologic_M_percentage_p_", p_value_str, ".png")
plot <- ggplot(pathologic_m_percentage, aes(x = WDR54, y = Percentage, fill = Pathologic_M)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#FF9999", "#FFCC99")) + 
  labs(title = "Percentage Stacked Bar Chart of Pathologic M by WDR54 Expression",
       x = "WDR54 Expression",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) +  
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 20, b = 15, l = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 1.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )
original_width <- 8  
new_width <- original_width * 0.618
ggsave(filename = output_path, plot = plot, width = new_width, height = 6, units = "in")
library(ggplot2)
library(dplyr)
library(reshape2)
file_path <- "C:/R/C_Project/F3/Clinical_landscape/Tumor_stage.csv"
tumor_stage_data <- read.csv(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
tumor_stage_long <- melt(tumor_stage_data, id.vars = "WDR54", variable.name = "Variable", value.name = "Tumor_stage")
tumor_stage_long <- tumor_stage_long[tumor_stage_long$Variable == "Tumor_stage",]
tumor_stage_long$Tumor_stage <- factor(tumor_stage_long$Tumor_stage, levels = c("1", "2", "3", "4"))
tumor_stage_percentage <- tumor_stage_long %>%
  group_by(WDR54, Tumor_stage) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

contingency_table <- table(tumor_stage_long$WDR54, tumor_stage_long$Tumor_stage)
chi_test <- chisq.test(contingency_table)
p_value <- chi_test$p.value
cat("P-value for the difference in Tumor_stage distribution between WDR54 groups:", p_value, "\n")
p_value_str <- sprintf("%.3e", p_value)
output_path <- paste0("C:/R/C_Project/F3/Clinical_landscape/WDR54_Tumor_stage_percentage_p_", p_value_str, ".png")
plot <- ggplot(tumor_stage_percentage, aes(x = WDR54, y = Percentage, fill = Tumor_stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#FF9999", "#FFCC99", "#99CCFF", "#66CC99")) + 
  labs(title = "Percentage Stacked Bar Chart of Tumor Stages by WDR54 Expression",
       x = "WDR54 Expression",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) +  
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 20, b = 15, l = 20),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, vjust = 1.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )
original_width <- 8 
new_width <- original_width * 0.618
ggsave(filename = output_path, plot = plot, width = new_width, height = 6, units = "in")
## KM curve
library(survival)
library(survminer)
setwd("C:/R/C_Project/F3/KM_CURVE")
clinical_data <- read.csv("CDMI_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$OS_time <- round(as.numeric(clinical_data$OS_time) * 30)
str(clinical_data)
write.csv(clinical_data, "CDMI_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = TRUE)
clinical_data <- read.csv("CDMI_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "CDMI_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "CDMI_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "CDMI_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE65858_PFS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(PFS_time, PFS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE65858_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE65858_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE65858_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE85446_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE85446_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE85446_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE85446_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE41613_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE41613_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE41613_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE41613_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE65858_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE65858_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE65858_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE65858_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE75538_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE75538_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE75538_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE75538_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("TCGA_HNSC_OS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(OS_time, OS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "TCGA_HNSC_OS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "TCGA_HNSC_OS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "TCGA_HNSC_OS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
library(survival)
library(survminer)
setwd("C:/R/C_Project/F3/KM_CURVE")
clinical_data <- read.csv("GSE75538_DFS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$DFS_time <- round(as.numeric(clinical_data$DFS_time) * 30)
str(clinical_data)
write.csv(clinical_data, "GSE75538_DFS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = TRUE)
clinical_data <- read.csv("GSE75538_DFS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(DFS_time, DFS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)
ggsave(filename = "GSE75538_DFS_survival_plot.png", plot = plot$plot, height = 6, width = 8)
ggsave(filename = "GSE75538_DFS_risk_table.png", plot = plot$table, height = 2, width = 8)
ggsave(filename = "GSE75538_DFS_ncensor_plot.png", plot = plot$ncensor.plot, height = 2.5, width = 8)
clinical_data <- read.csv("GSE65858_PFS_FROM_CDMI_CLINICAL_ANNOTATION.csv", row.names = 1)
str(clinical_data)
clinical_data$gene_group <- ifelse(clinical_data$WDR54 > median(clinical_data$WDR54, na.rm = TRUE), "High", "Low")
fit <- survfit(Surv(PFS_time, PFS) ~ gene_group, data = clinical_data)
plot <- ggsurvplot(
    fit,
    data = clinical_data,
    title = "Survival curves based on WDR54 expression", 
    subtitle = "Based on Kaplan-Meier estimates",
    caption = "created with survminer",
    pval = TRUE, 
    conf.int = TRUE, 
    font.title = c(16, "bold", "darkblue"),
    font.subtitle = c(15, "bold.italic", "purple"),
    font.caption = c(14, "plain", "orange"),
    font.x = c(14, "bold.italic", "red"),
    font.y = c(14, "bold.italic", "darkred"),
    font.tickslab = c(12, "plain", "darkgreen"),
    risk.table = TRUE,
    risk.table.title = "Note the risk set sizes",
    risk.table.subtitle = "and remember about censoring.",
    risk.table.caption = "source code: website.com",
    risk.table.height = 0.2,
    ncensor.plot = TRUE,
    ncensor.plot.title = "Number of censorings",
    ncensor.plot.subtitle = "over the time.",
    ncensor.plot.caption = "data available at data.com",
    ncensor.plot.height = 0.25,
    height = 4  
)
print(plot)

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
rm(list = ls())
options(stringsAsFactors = FALSE)

file_path <- "C:/R/C_Project/F5/WDR54_compare/GSE65858_Subtype.csv"
plotdata <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
head(plotdata)

plotdata <- plotdata %>%
  dplyr::rename(Sample = X, WDR54 = WDR54, Group = Tissue_type) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Mesenchymal", "Atypical", "Basal", "Classical")))

comparisons <- list(
  c("Mesenchymal", "Atypical"),
  c("Mesenchymal", "Basal"),
  c("Mesenchymal", "Classical")
)
p_values <- compare_means(WDR54 ~ Group, data = plotdata, comparisons = comparisons, method = "t.test")

p_values_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_Subtype_P_values.csv"
write.csv(p_values, file = p_values_path, row.names = FALSE)

pl <- ggplot(data = plotdata, aes(x = Group, y = WDR54)) +
  geom_violin(aes(fill = Group), alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.1, notch = TRUE) +
  geom_quasirandom(aes(color = Group), shape = 16, size = 1, width = 0.1) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  labs(x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#AF0F11", "#3372A6", "#367B34", "#D95F02")) +
  scale_fill_manual(values = c("#AF0F11", "#3372A6", "#367B34", "#D95F02")) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_Subtype_Plot.png"
ggsave(filename = output_path, plot = pl, width = 8, height = 6, units = "in", dpi = 300)

install.packages("pROC")
install.packages("ggplot2")
install.packages("dplyr")

library(pROC)
library(ggplot2)
library(dplyr)

file_path <- "C:/R/C_Project/F5/WDR54_compare/GSE65858_Subtype.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

data <- data %>%
  dplyr::mutate(Group = ifelse(Tissue_type == "Mesenchymal", "Mesenchymal", "Non-Mesenchymal"))
data$Group <- factor(data$Group, levels = c("Mesenchymal", "Non-Mesenchymal"))

roc_data <- roc(data$Group, data$WDR54)

png(
  filename = "C:/R/C_Project/F5/WDR54_compare/ROC_Curve_Plot.png",
  width    = 2400,
  height   = 1800,
  res      = 300
)
plot(
  roc_data,
  col          = rgb(196/255, 112/255, 112/255),
  legacy.axes  = TRUE,
  main         = "ROC Curve for WDR54 Expression",
  xlab         = "1 - Specificity",
  ylab         = "Sensitivity",
  lwd          = 2,
  xlim         = c(0, 1),
  ylim         = c(0, 1),
  print.auc    = FALSE,
  print.thres  = FALSE
)
plot(
  roc_data,
  col               = rgb(196/255, 112/255, 112/255),
  legacy.axes       = TRUE,
  auc.polygon       = TRUE,
  grid              = c(0.1, 0.2),
  grid.col          = c("green", "red"),
  max.auc.polygon   = TRUE,
  auc.polygon.col   = rgb(219/255, 184/255, 178/255),
  print.thres       = FALSE,
  xlim              = c(0, 1),
  ylim              = c(0, 1)
)
dev.off()

auc_value <- round(auc(roc_data), 3)
ci_value  <- round(ci(roc_data), 3)

results <- data.frame(
  Measure = c("AUC", "95% CI Lower", "95% CI Upper"),
  Value   = c(auc_value, ci_value[1], ci_value[2])
)
write.csv(results, file = "C:/R/C_Project/F5/WDR54_compare/ROC_Curve_Results.csv", row.names = FALSE)
cat("ROC AUC:", auc_value, "95% CI:", ci_value, "\n")

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
rm(list = ls())
options(stringsAsFactors = FALSE)

file_path <- "C:/R/C_Project/F5/WDR54_compare/GSE201777_Primary_VS_Metastatic.csv"
plotdata <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
head(plotdata)

plotdata <- plotdata %>%
  dplyr::rename(Sample = X, WDR54 = WDR54, Group = Tissue_type) %>%
  dplyr::mutate(Group = factor(Group))

comparisons <- list(c("Primary", "Metastatic"))
p_values   <- compare_means(WDR54 ~ Group, data = plotdata, comparisons = comparisons, method = "t.test")
custom_colors <- c("#AF0F11", "#3372A6")
p_values_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_GSE201777_Primary_VS_Metastatic_P_values.csv"
write.csv(p_values, file = p_values_path, row.names = FALSE)

pl <- ggplot(data = plotdata, aes(x = Group, y = WDR54)) +
  geom_violin(aes(fill = Group), alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.1, notch = FALSE) +
  geom_quasirandom(aes(color = Group), shape = 16, size = 2, width = 0.1) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  labs(x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_GSE201777_Primary_VS_Metastatic_Plot.png"
ggsave(filename = output_path, plot = pl, width = 5, height = 6, units = "in", dpi = 300)

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
rm(list = ls())
options(stringsAsFactors = FALSE)

file_path <- "C:/R/C_Project/F5/WDR54_compare/TCGA_Invasion.csv"
plotdata <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
plotdata <- plotdata %>%
  dplyr::mutate(Merged = factor(Merged, levels = c("YES", "NO")))

pl <- ggplot(data = plotdata, aes(x = Merged, y = WDR54)) +
  geom_violin(aes(fill = Merged), alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.1, notch = TRUE) +
  geom_quasirandom(aes(color = Merged), shape = 16, size = 1, width = 0.1) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(x = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#AF0F11", "#3372A6")) +
  scale_fill_manual(values = c("#AF0F11", "#3372A6")) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_Invasion_Plot.png"
ggsave(filename = output_path, plot = pl, width = 5, height = 6, units = "in", dpi = 300)

rm(list = ls())
options(stringsAsFactors = FALSE)

input_file <- "C:/R/C_Project/F5/WDR54_compare/TCGA_PN.csv"
plotdata <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
plotdata$pathologic_N <- ifelse(plotdata$pathologic_N == "N0", "NO", "YES")
plotdata$pathologic_N <- factor(plotdata$pathologic_N, levels = c("YES", "NO"))

pl <- ggplot(data = plotdata, aes(x = pathologic_N, y = WDR54)) +
  geom_violin(aes(fill = pathologic_N), alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.1, notch = TRUE) +
  geom_quasirandom(aes(color = pathologic_N), shape = 16, size = 1, width = 0.1) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(x = "Pathologic N (Binary)", y = "WDR54 Expression") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#AF0F11", "#3372A6")) +
  scale_fill_manual(values = c("#AF0F11", "#3372A6")) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_PN_Binary_Plot.png"
ggsave(filename = output_path, plot = pl, width = 5, height = 6, units = "in", dpi = 300)

rm(list = ls())
options(stringsAsFactors = FALSE)

input_file <- "C:/R/C_Project/F5/WDR54_compare/TCGA_Grade.csv"
plotdata <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
plotdata$neoplasm_histologic_grade <- factor(plotdata$neoplasm_histologic_grade, levels = c("G1", "G2", "G3"))

pl <- ggplot(data = plotdata, aes(x = neoplasm_histologic_grade, y = WDR54)) +
  geom_violin(aes(fill = neoplasm_histologic_grade), alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.1, notch = TRUE) +
  geom_quasirandom(aes(color = neoplasm_histologic_grade), shape = 16, size = 1, width = 0.1) +
  labs(x = "Neoplasm Histologic Grade (G1-G3)", y = "WDR54 Expression") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(values = c("#AF0F11", "#3372A6", "#367B34")) +
  scale_fill_manual(values = c("#AF0F11", "#3372A6", "#367B34")) +
  guides(color = "none", shape = "none", fill = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text  = element_text(size = 10, color = "black"),
    text       = element_text(size = 9, color = "black")
  )

output_path <- "C:/R/C_Project/F5/WDR54_compare/WDR54_Grade_Plot.png"
ggsave(filename = output_path, plot = pl, width = 6.5, height = 6, units = "in", dpi = 300)

library(pROC)
library(ggplot2)
library(dplyr)

file_path <- "C:/R/C_Project/F5/WDR54_compare/TCGA_Grade.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
data <- data %>%
  dplyr::mutate(Group = ifelse(neoplasm_histologic_grade == "G1", "Well differentiated", "Poorly differentiated"))
data$Group <- factor(data$Group, levels = c("Well differentiated", "Poorly differentiated"))

roc_data <- roc(data$Group, data$WDR54)

png(filename = "C:/R/C_Project/F5/WDR54_compare/ROC_Curve_Grade_Plot.png", width = 2400, height = 1800, res = 300)
plot(roc_data, col = "#5CB3BA", legacy.axes = TRUE, main = "ROC Curve for WDR54 Expression (Well vs Poorly differentiated)", xlab = "1 - Specificity", ylab = "Sensitivity", lwd = 2, xlim = c(0, 1), ylim = c(0, 1), print.auc = FALSE, print.thres = FALSE)
plot(roc_data, col = "#5CB3BA", legacy.axes = TRUE, auc.polygon = TRUE, grid = c(0.1, 0.2), grid.col = c("green", "red"), max.auc.polygon = TRUE, auc.polygon.col = "#5CB3BA", print.thres = FALSE, xlim = c(0, 1), ylim = c(0, 1))
dev.off()

auc_value <- round(auc(roc_data), 3)
ci_value  <- round(ci(roc_data), 3)
results <- data.frame(Measure = c("AUC", "95% CI Lower", "95% CI Upper"), Value = c(auc_value, ci_value[1], ci_value[2]))
write.csv(results, file = "C:/R/C_Project/F5/WDR54_compare/ROC_Curve_Grade_Results.csv", row.names = FALSE)
cat("ROC AUC:", auc_value, "95% CI:", ci_value, "\n")

## survival analysis
library(survival)
library(forestplot)
library(dplyr)

file_path <- "C:/R/C_Project/SF5/TCGA_COX.csv"
data <- read.csv(file_path)

cox_model <- coxph(Surv(OS_time, OS) ~ WDR54 + Age + Gender + Tumor_Grade, data = data)
summary_cox <- summary(cox_model)

result <- data.frame(
  HR = round(summary_cox$conf.int[, "exp(coef)"], 2),
  Lower_CI = round(summary_cox$conf.int[, "lower .95"], 2),
  Upper_CI = round(summary_cox$conf.int[, "upper .95"], 2),
  P_value = formatC(summary_cox$coefficients[, "Pr(>|z|)"], format = "f", digits = 3),
  Characteristics = rownames(summary_cox$conf.int)
)

output_file <- "C:/R/C_Project/SF5/COX_forest_plot.png"

png(filename = output_file, width = 800, height = 600)
forestplot(
  result[, c("Characteristics", "HR", "P_value")],
  mean  = as.numeric(result$HR),
  lower = as.numeric(result$Lower_CI),
  upper = as.numeric(result$Upper_CI),
  zero      = 1,
  boxsize   = 0.2,
  graph.pos = 2,
  title     = "Single Variable Cox Analysis for Multiple Variables",
  col       = fpColors(box = "#021eaa", lines = "#021eaa", zero = "black"),
  txt_gp    = fpTxtGp(
    label   = gpar(cex = 1.5),
    ticks   = gpar(cex = 1.5),
    xlab    = gpar(cex = 2),
    title   = gpar(cex = 2.5),
    summary = gpar(cex = 1.5)
  ),
  fn.ci_norm = "fpDrawDiamondCI",
  txt_gp_labels = list(
    label = list(
      gpar(col = rep("black", nrow(result)))
    )
  )
)
dev.off()

cat("Cox回归分析和森林图已成功生成，并保存至:", output_file, "\n")

library(survival)
library(forestplot)
library(tableone)

file_path <- "C:/R/C_Project/SF5/TCGA_COX.csv"
data <- read.csv(file_path)
str(data)

result_list <- list()
variables <- c("WDR54", "Age", "Gender", "Tumor_Grade")

for (var in variables) {
  cox_model   <- coxph(as.formula(paste("Surv(OS_time, OS) ~", var)), data = data)
  cox_summary <- summary(cox_model)
  hr          <- round(cox_summary$conf.int[1, 1], 2)
  lower_ci    <- round(cox_summary$conf.int[1, 3], 2)
  upper_ci    <- round(cox_summary$conf.int[1, 4], 2)
  p_value     <- round(cox_summary$coefficients[1, 5], 3)
  result      <- data.frame(
    Characteristics = var,
    HR             = hr,
    Lower_CI       = lower_ci,
    Upper_CI       = upper_ci,
    P_value        = p_value
  )
  result_list[[var]] <- result
}

result <- do.call(rbind, result_list)
print(result)

output_file <- "C:/R/C_Project/SF5/Single_Variable_Cox_Forest_Plot.png"
png(filename = output_file, width = 800, height = 600)
forestplot(
  result[, c("Characteristics", "HR", "P_value")],
  mean       = result$HR,
  lower      = result$Lower_CI,
  upper      = result$Upper_CI,
  zero       = 1,
  boxsize    = 0.2,
  graph.pos  = 2,
  title      = "Single Variable Cox Analysis for Multiple Variables",
  col        = fpColors(box = "#021eaa", lines = "#021eaa", zero = "black"),
  txt_gp     = fpTxtGp(
                  label   = gpar(cex = 1.5),
                  ticks   = gpar(cex = 1.5),
                  xlab    = gpar(cex = 2),
                  title   = gpar(cex = 2.5),
                  summary = gpar(cex = 1.5)
                ),
  fn.ci_norm     = "fpDrawDiamondCI",
  txt_gp_labels  = list(
                     label = list(
                       gpar(col = rep("black", nrow(result)))
                     )
                   )
)

dev.off()

rm(list = ls())
library(survival)
library(survminer)
library(rms)

file_path <- "C:/R/C_Project/SF5/TCGA_COX.csv"
data <- read.csv(file_path)

data$Gender      <- factor(data$Gender, levels = c(1, 2), labels = c("Male", "Female"))
data$Tumor_Grade <- factor(data$Tumor_Grade,
                           levels = c(1, 2, 3, 4),
                           labels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"))

dd <- datadist(data)
options(datadist = "dd")

f2 <- psm(Surv(OS_time, OS) ~ WDR54 + Age + Gender + Tumor_Grade,
          data = data, dist = "lognormal")
f2 <- psm(Surv(OS_time, OS) ~ WDR54 + Age + Gender + Tumor_Grade,
          data = data, x = TRUE, y = TRUE, dist = "lognormal")

cal1 <- calibrate(f2, cmethod = "KM", method = "boot", u = 365,  m = 76, B = 1000)
cal3 <- calibrate(f2, cmethod = "KM", method = "boot", u = 1095, m = 76, B = 1000)
cal5 <- calibrate(f2, cmethod = "KM", method = "boot", u = 1825, m = 76, B = 1000)

png(filename = "C:/R/C_Project/SF5/calibration_curve_1_3_5_years.png",
    width    = 800,
    height   = 600)

plot(cal1, lwd = 2, lty = 0, errbar.col = "#2166AC",
     bty      = "l",
     xlim     = c(0, 1), ylim = c(0, 1),
     xlab     = "Nomogram-predicted OS (%)",
     ylab     = "Observed OS (%)",
     col      = "#2166AC",
     cex.lab  = 1.2,
     cex.axis = 1,
     cex.main = 1.2,
     cex.sub  = 0.6)

lines(cal1[, c("mean.predicted", "KM")], type = "b", lwd = 1, col = "#2166AC", pch = 16)
lines(cal3[, c("mean.predicted", "KM")], type = "b", lwd = 2, col = "#B2182B", pch = 16)
lines(cal5[, c("mean.predicted", "KM")], type = "b", lwd = 2, col = "#F4A582", pch = 16)

abline(0, 1, lwd = 2, lty = 3, col = "#224444")

legend("topleft",
       legend = c("1-year", "3-year", "5-year"),
       col    = c("#2166AC", "#B2182B", "#F4A582"),
       lwd    = 2,
       cex    = 1.2,
       bty    = "n")

dev.off()
