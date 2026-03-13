
## setup
library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(Matrix)
library(patchwork)
library(Nebulosa)
library(ggrepel)
library(cowplot)
library(circlize)
library(fgsea)
library(ggrastr)
library(ggpubr)
library(tidyr)



##############################  import  ############################################### 
setwd("/Users/birzh586/Desktop/skin-horse/analysis")
merge <- readRDS("/Users/birzh586/Desktop/skin-horse/data/merge.rds")
merge
# An object of class Seurat 
# 30152 features across 100434 samples within 1 assay 
# Active assay: RNA (30152 features, 0 variable features)
# 1 layer present: counts

head(merge)
# orig.ident nCount_RNA nFeature_RNA sample
# C241029001a_TGCCTGATC_AACGCTAGT_AACAAGTGG  TGCCTGATC       3602         1376   C_1a
# C241029001a_ACATCGGAC_AACGTCCAA_AACAAGTGG  ACATCGGAC        691          419   C_1a
# C241029001a_CCACACATT_AACGTCCAA_AACAAGTGG  CCACACATT        862          495   C_1a

table(merge$sample)
# C_1a  C_1b  C_2a  C_2b 
# 11680 31878 24989 31887 
############################## QC ################################
dir.create("1.quality_control")
setwd("/Users/birzh586/Desktop/skin-horse/analysis/1.quality_control")

# DS_masto_integrated_PCA.rds
DefaultAssay(merge) <- "RNA"

# store mitochondrial percentage in object meta data
merge
merge <- PercentageFeatureSet(merge, pattern = "^MT.", col.name = "percent.mt")

# Visualize QC metrics as a violin plot
pdf("QC-plot.pdf", width = 18, height = 3)
VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, group.by = "sample", raster=T)
dev.off()

# filter
# merge <- subset(merge, subset = nFeature_RNA > 200 & nCount_RNA < 10000 & percent.mt < 60)

## metadata
metadata <- merge@meta.data 
head(metadata)

summary_df <- metadata %>%
  group_by(sample) %>%
  summarise(
    nCount_q5 = quantile(nCount_RNA, 0.05, na.rm = TRUE),
    nCount_q95 = quantile(nCount_RNA, 0.95, na.rm = TRUE),
    
    nFeature_q5 = quantile(nFeature_RNA, 0.05, na.rm = TRUE),
    nFeature_q95 = quantile(nFeature_RNA, 0.95, na.rm = TRUE),
    
    percent_mt_q5 = quantile(percent.mt, 0.05, na.rm = TRUE),
    percent_mt_q95 = quantile(percent.mt, 0.95, na.rm = TRUE)
  )
write.csv(summary_df, "summary_df_raw.csv",row.names = F)

# summary_df <- metadata %>%
#   group_by(sample) %>%
#   summarise(
#     nCount_min = min(nCount_RNA, na.rm = TRUE),
#     nCount_median = median(nCount_RNA, na.rm = TRUE),
#     nCount_max = max(nCount_RNA, na.rm = TRUE),
#     nCount_q5 = quantile(nCount_RNA, 0.05, na.rm = TRUE),
#     nCount_q95 = quantile(nCount_RNA, 0.95, na.rm = TRUE),
#     
#     nFeature_min = min(nFeature_RNA, na.rm = TRUE),
#     nFeature_median = median(nFeature_RNA, na.rm = TRUE),
#     nFeature_max = max(nFeature_RNA, na.rm = TRUE),
#     nFeature_q5 = quantile(nFeature_RNA, 0.05, na.rm = TRUE),
#     nFeature_q95 = quantile(nFeature_RNA, 0.95, na.rm = TRUE),
#     
#     percent_mt_min = min(percent.mt, na.rm = TRUE),
#     percent_mt_median = median(percent.mt, na.rm = TRUE),
#     percent_mt_max = max(percent.mt, na.rm = TRUE),
#     percent_mt_q5 = quantile(percent.mt, 0.05, na.rm = TRUE),
#     percent_mt_q95 = quantile(percent.mt, 0.95, na.rm = TRUE)
#   )


# # Define quantile probabilities (0, 0.1, ..., 1)
# probs <- seq(0, 1, 0.1)
# 
# # Function to compute quantiles per sample for a given metric, then pivot wide
# compute_quantiles_wide <- function(data, variable) {
#   data %>%
#     group_by(sample) %>%
#     summarise(
#       Quantiles = list(quantile(.data[[variable]], probs = probs, na.rm = TRUE))
#     ) %>%
#     ungroup() %>%
#     mutate(Quantile = list(probs)) %>%
#     unnest(cols = c(Quantile, Quantiles)) %>%
#     pivot_wider(
#       names_from = sample,
#       values_from = Quantiles
#     )
# }
# 
# # Compute quantiles for each metric (each will have Quantile + one column per sample)
# qc_nCount_wide <- compute_quantiles_wide(metadata, "nCount_RNA")
# qc_nFeature_wide <- compute_quantiles_wide(metadata, "nFeature_RNA")
# qc_percent_mt_wide <- compute_quantiles_wide(metadata, "percent.mt")
# 
# # View results
# qc_nCount_wide
# write.csv(qc_nCount_wide, "qc_nCount_wide.csv", row.names = F)
# qc_nFeature_wide
# write.csv(qc_nFeature_wide, "qc_nFeature_wide.csv",row.names = F)
# qc_percent_mt_wide
# qc_percent_mt_wide[,-1] <- round(qc_percent_mt_wide[,-1],4)
# write.csv(qc_percent_mt_wide, "qc_percent_mt_wide.csv",row.names = F)
# 


pdf("nCount_RNA_raw.pdf", width = 3, height = 3)
p1<- ggplot(metadata, aes(x = sample, y = nCount_RNA, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 0.5
  )+
  # 箱型图
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    color = "grey50",
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "grey50",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  scale_color_manual(values= scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  labs(
    x = "",
    y = "",
    title = "nCount_RNA"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p1
dev.off()





head(metadata)
pdf("nFeature_RNA_raw.pdf", width = 3, height = 3)
p2 <- ggplot(metadata, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 0.5
  )+
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    color = "grey50",
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "grey50",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  scale_color_manual(values= scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  labs(
    x = "",
    y = "",
    title = "nFeature_RNA"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p2
dev.off()





head(metadata)
pdf("percent.mt_raw.pdf", width = 3, height = 3)
p3<-ggplot(metadata, aes(x = sample, y = percent.mt, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 0.5
  )+
  # 箱型图
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    color = "grey50",
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "grey50",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  scale_color_manual(values= scales::alpha(c("bisque","darksalmon","#F67280","#C06C84"), 0.05))+
  labs(
    x = "",
    y = "",
    title = "Percent of MT"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p3
dev.off()



pdf("Quality Control_raw.pdf", width = 9, height = 3)
# 组合图
combined <- ggarrange(
  p1, p2, p3,
  widths = c(1, 1., 0.95),
  ncol = 3, nrow = 1
)
combined
dev.off()







# 
# # 创建一个空向量存储保留的细胞
# cells_to_keep <- c()
# 
# # 遍历每个样本
# for(s in unique(merge$sample)) {
#   
#   cat("Processing sample:", s, "\n")
#   
#   # 取该样本细胞
#   # 使用元数据列筛选细胞
#   sample_cells <- rownames(merge@meta.data)[merge$sample == s]
#   sample_data <- merge@meta.data[sample_cells, ]
#   
#   # 使用 IQR 计算上下界
#   get_bounds <- function(x, lower_quant = 0.05, upper_quant = 0.95, factor = 1.5) {
#     q_low <- quantile(x, lower_quant)
#     q_high <- quantile(x, upper_quant)
#     iqr <- q_high - q_low
#     lower <- q_low - factor * iqr
#     upper <- q_high + factor * iqr
#     return(c(lower, upper))
#   }
#   
#   # nFeature_RNA
#   bounds_nFeature <- get_bounds(sample_data$nFeature_RNA)
#   # nCount_RNA
#   bounds_nCount <- get_bounds(sample_data$nCount_RNA)
#   # percent.mt (只过滤上界)
#   upper_mt <- quantile(sample_data$percent.mt, 0.95) + 1.5 * IQR(sample_data$percent.mt)
#   # upper_mt <- 2
#   
#   # 筛选符合条件的细胞
#   keep <- rownames(sample_data)[
#     sample_data$nFeature_RNA > bounds_nFeature[1] &
#       sample_data$nFeature_RNA < bounds_nFeature[2] &
#       sample_data$nCount_RNA > bounds_nCount[1] &
#       sample_data$nCount_RNA < bounds_nCount[2] &
#       sample_data$percent.mt < upper_mt
#   ]
#   
#   cells_to_keep <- c(cells_to_keep, keep)
# }



# ===================================
# 方法1: 混合策略（推荐）⭐
# ===================================
adaptive_qc_filter <- function(seurat_obj, sample_col = "sample") {
  
  # 获取所有样本
  samples <- unique(seurat_obj@meta.data[[sample_col]])
  
  # 存储每个样本的cutoff
  cutoffs_list <- list()
  filtered_cells <- c()
  
  for (sample in samples) {
    # 提取该样本的metadata
    sample_meta <- seurat_obj@meta.data[seurat_obj@meta.data[[sample_col]] == sample, ]
    
    # 1. nFeature_RNA: 
    nFeature_lower <- 200  # Biological lower limit
    nFeature_upper <- min(
      quantile(sample_meta$nFeature_RNA, 0.95),
      6000  # Biological upper limit
    )
    
    # 2. nCount_RNA: Median Absolute Deviation
    median_count <- median(sample_meta$nCount_RNA)
    mad_count <- mad(sample_meta$nCount_RNA)
    nCount_lower <- max(
      median_count - 3 * mad_count,
      500
    )
    nCount_upper <- median_count + 5 * mad_count
    
    # 3. percent.mt
    mt_median <- median(sample_meta$percent.mt)
    mt_mad <- mad(sample_meta$percent.mt)
    
    if (mt_median < 1) {
      mt_cutoff <- 5.0  # Standard threshold
    } else {
      mt_cutoff <- min(
        mt_median + 3 * mt_mad,
        20.0
      )
    }
    
    # 应用过滤
    pass_qc <- (
      sample_meta$nFeature_RNA > nFeature_lower &
        sample_meta$nFeature_RNA < nFeature_upper &
        sample_meta$nCount_RNA > nCount_lower &
        sample_meta$nCount_RNA < nCount_upper &
        sample_meta$percent.mt < mt_cutoff
    )
    
    filtered_cells <- c(filtered_cells, rownames(sample_meta)[pass_qc])
    
    # 保存cutoff信息
    cutoffs_list[[sample]] <- list(
      nFeature = c(nFeature_lower, nFeature_upper),
      nCount = c(nCount_lower, nCount_upper),
      mt = mt_cutoff,
      before = nrow(sample_meta),
      after = sum(pass_qc),
      retained_pct = round(sum(pass_qc) / nrow(sample_meta) * 100, 2)
    )
    
    # 打印信息
    cat("\n", sample, ":\n", sep = "")
    cat("  nFeature: ", round(nFeature_lower), " - ", round(nFeature_upper), "\n", sep = "")
    cat("  nCount: ", round(nCount_lower), " - ", round(nCount_upper), "\n", sep = "")
    cat("  percent.mt: < ", round(mt_cutoff, 2), "%\n", sep = "")
    cat("  Before: ", nrow(sample_meta), " cells\n", sep = "")
    cat("  After: ", sum(pass_qc), " cells (", round(sum(pass_qc)/nrow(sample_meta)*100, 1), "%)\n", sep = "")
  }
  
  # 过滤Seurat对象
  merge_filtered <- subset(seurat_obj, cells = filtered_cells)
  
  cat("\n=== Summary ===\n")
  cat("Total before: ", ncol(seurat_obj), "\n")
  cat("Total after: ", ncol(merge_filtered), "\n")
  cat("Retained: ", round(ncol(merge_filtered)/ncol(seurat_obj)*100, 2), "%\n")
  
  return(list(
    seurat = merge_filtered,
    cutoffs = cutoffs_list
  ))
}

# use adaptive_qc_filter() function
result <- adaptive_qc_filter(merge, sample_col = "sample")
merge_filtered <- result$seurat
cutoffs_info <- result$cutoffs


print(cutoffs_info)


merge_filtered
# An object of class Seurat 
# 30152 features across 93002 samples within 1 assay 
# Active assay: RNA (30152 features, 0 variable features)
# 1 layer present: counts


table(merge_filtered@meta.data$sample)
# C_1a C_1b C_2a C_2b 
#    10263       29539       23508       29692 
saveRDS(merge_filtered, file = "/Users/birzh586/Desktop/skin-horse/data/merge_filtered.rds")


## metadat
metadata_filter <- merge_filtered@meta.data 
head(metadata_filter)
table(metadata_filter$sample)



summary_df_filter <- metadata_filter %>%
  group_by(sample) %>%
  summarise(
    nCount_q5 = quantile(nCount_RNA, 0.05, na.rm = TRUE),
    nCount_q95 = quantile(nCount_RNA, 0.95, na.rm = TRUE),
    
    nFeature_q5 = quantile(nFeature_RNA, 0.05, na.rm = TRUE),
    nFeature_q95 = quantile(nFeature_RNA, 0.95, na.rm = TRUE),
    
    percent_mt_q5 = quantile(percent.mt, 0.05, na.rm = TRUE),
    percent_mt_q95 = quantile(percent.mt, 0.95, na.rm = TRUE)
  )

write.csv(summary_df_filter,"summary_df_filter.csv")









pdf("nCount_RNA_filter.pdf", width = 3, height = 3)
p4<- ggplot(metadata_filter, aes(x = sample, y = nCount_RNA, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 1
  )+
  # 箱型图
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "black",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  scale_color_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  labs(
    x = "",
    y = "",
    title = "nCount_RNA"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p4
dev.off()





head(metadata_filter)
pdf("nFeature_RNA_filter.pdf", width = 3, height = 3)
p5 <- ggplot(metadata_filter, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 1
  )+
  # 箱型图
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "black",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  scale_color_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  labs(
    x = "",
    y = "",
    title = "nFeature_RNA"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p5
dev.off()





head(metadata_filter)
pdf("percent.mt_filter.pdf", width = 3, height = 3)
p6<-ggplot(metadata_filter, aes(x = sample, y = percent.mt, fill = sample)) +
  # 每个样本点
  geom_jitter_rast(
    aes(color = sample),
    position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8),
    size = 0.1,
    alpha = 1
  )+
  # 箱型图
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 1,
    outlier.shape = NA  # 不显示箱型图自身的离群点
  ) +
  # 横线：中位数
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.2,
    position = position_dodge(width = 0.55),
    color = "black",
    size = 0.1
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  scale_color_manual(values=c("bisque","darksalmon","#F67280","#C06C84"))+
  labs(
    x = "",
    y = "",
    title = "Percent of MT"
  ) +
  theme_bw()+  
  theme(panel.grid=element_blank(), 
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill= "white"),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, color="black",hjust = 1, angle=60),
        axis.text.y = element_text(size = 11, color="black"),
        axis.ticks.x = element_blank(), 
        legend.title =  element_blank(), 
        legend.text = element_text(size = 13),
        legend.position="none")
p6
dev.off()



pdf("Quality Control_filter.pdf", width = 9, height = 3)
# 组合图
combined <- ggarrange(
  p4, p5, p6,
  widths = c(1, 1, 0.95),
  ncol = 3, nrow = 1
)
combined
dev.off()







color <- c("bisque","darksalmon","#F67280","#C06C84")
##################################### old ###################################
table(merge$sample)
# C_1a C_1b C_2a C_2b 
# 11680       31878       24989       31887 

11680    +   31878    +   24989   +    31887 
# 100434
old <- data.frame(
  samples = c("C_1a", "C_1b", "C_2a", "C_2b"),
  Freq = c(11680, 31878, 24989, 31887)
)

##fraction
old$fraction <- old$Freq/sum(old$Freq)
old$fraction <- as.numeric(old$fraction)
# Compute the cumulative percentages (top of each rectangle)
old$ymax <- cumsum(old$fraction)
# Compute the bottom of each rectangle
old$ymin <- c(0, head(old$ymax, n=-1))

## category to factor
level <- old$samples
old$samples <- factor(old$samples , levels=level)
old$samples



head(old)
# 计算 label 位置
old <- old %>%
  #arrange(desc(label)) %>%
  mutate(
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2,
    percent_label = paste0(round(fraction * 100,2), "%")
  )


head(old)
pdf("sample_size_raw.pdf", width = 5, height = 3)
p7 <- ggplot(old, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = samples)) +
  geom_rect(color = "grey50", size = 0.5) +
  coord_polar(theta = "y") +
  theme_void() +
  xlim(c(0.1, 5)) +  # 给标签留出足够空间
  geom_text(
    aes(x = 3, y = label_pos, label = Freq),  # x 选在条形中间
    size = 4
  )+
  ggtitle("Raw = 100,434 cells")+
  scale_fill_manual(values = color) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12, hjust = 0.5, vjust=-10),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank() ,   
    legend.background = element_blank(), 
    legend.key        = element_blank()
  ) +
  guides(fill = guide_legend(color = "none"))
p7
dev.off()










##################################### new ###################################
table(merge_filtered$sample)
# C_1a C_1b C_2a C_2b 
#  10263       29539       23508       29692 

10263      + 29539    +   23508    +   29692 
# 93002

new <- data.frame(
  samples = c("C_1a", "C_1b", "C_2a", "C_2b"),
  Freq = c(10263,29539,23508,29692 )
)

##fraction
new$fraction <- new$Freq/sum(new$Freq)
new$fraction <- as.numeric(new$fraction)
# Compute the cumulative percentages (top of each rectangle)
new$ymax <- cumsum(new$fraction)
# Compute the bottom of each rectangle
new$ymin <- c(0, head(new$ymax, n=-1))

## category to factor
level <- new$samples
new$samples <- factor(new$samples , levels=level)
new$samples



head(new)
# 计算 label 位置
new <- new %>%
  #arrange(desc(label)) %>%
  mutate(
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2,
    percent_label = paste0(round(fraction * 100,2), "%")
  )


head(new)
pdf("sample_size_filter.pdf", width = 5, height = 3)
p8 <- ggplot(new, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = samples)) +
  geom_rect(color = "black", size = 0.5) +
  coord_polar(theta = "y") +
  theme_void() +
  xlim(c(0.1, 5)) +  # 给标签留出足够空间
  geom_text(
    aes(x = 3, y = label_pos, label = Freq),  # x 选在条形中间
    size = 4
  )+
  ggtitle("Filter = 93,002 cells")+
  scale_fill_manual(values = color) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12, hjust = 0.5, vjust=-10),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank() ,   
    legend.background = element_blank(), 
    legend.key        = element_blank()
  ) +
  guides(fill = guide_legend(color = "none"))
p8
dev.off()




##################################### final ###################################
pdf("Quality Control_final.pdf", width = 12, height = 6)
# 组合图
combined <- ggarrange(
  p7,p1,p2,p3,
  p8,p4,p5,p6,
  widths = c(1.1, 1, 1, 0.95),
  ncol = 4, nrow = 2
)
combined
dev.off()


