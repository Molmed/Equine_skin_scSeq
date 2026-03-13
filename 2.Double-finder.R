setwd("/Users/birzh586/Desktop/skin-horse/analysis")
dir.create("2.doubletfinder")
setwd("/Users/birzh586/Desktop/skin-horse/analysis/2.doubletfinder")


seu <- readRDS("/Users/birzh586/Desktop/skin-horse/data/merge_filtered.rds")
seu
# An object of class Seurat 
# 30152 features across 93002 samples within 1 assay 
# Active assay: RNA (30152 features, 0 variable features)
# 1 layer present: counts

library(Seurat)
library(DoubletFinder)
library(BiocNeighbors)
library(BiocParallel)
library(dplyr)
library(tibble)

# ============================================================
# 完整修复版 run_doubletfinder_custom
# ============================================================
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL) {
  
  sample_name <- unique(seu_sample_subset[["sample"]])
  print(paste0("Processing sample: ", sample_name, "..."))
  
  multiplet_rates_10x <- data.frame(
    Multiplet_rate  = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    Recovered_cells = c(500,   1000,  2000,  3000,  4000,  5000,  6000,  7000,  8000,  9000,  10000)
  )
  
  
  # ============================================================
  # 1. 估算 doublet 率
  # ============================================================
  if (is.null(multiplet_rate)) {
    
    nCells <- nrow(seu_sample_subset@meta.data)
    
    # 替换查表部分，改为线性插值+合理上限
    if (nCells <= 500) {
      multiplet_rate <- 0.004
    } else {
      lm_model <- lm(Multiplet_rate ~ Recovered_cells, data = multiplet_rates_10x)
      multiplet_rate <- min(
        as.numeric(predict(lm_model, newdata = data.frame(Recovered_cells = nCells))),
        0.10  # 上限改为10%，更保守
      )
    }
    cat("细胞数:", nCells, "| doublet率:", round(multiplet_rate * 100, 1), "%\n")
  }
  
  # ============================================================
  # 2. 预处理
  # ============================================================
  sample <- NormalizeData(seu_sample_subset, verbose = FALSE)
  sample <- FindVariableFeatures(sample, verbose = FALSE)
  sample <- ScaleData(sample, verbose = FALSE)
  sample <- RunPCA(sample, verbose = FALSE)
  
  # 自动选最优 PC 数
  stdv         <- sample[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative   <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[1:(length(percent_stdv)-1)] -
                       percent_stdv[2:length(percent_stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)
  min_pc <- max(min_pc, 5)   # 至少用5个PC
  min_pc <- min(min_pc, 30)  # 最多用30个PC
  cat("使用 PC 数:", min_pc, "\n")
  
  sample <- RunUMAP(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindNeighbors(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindClusters(sample, resolution = 0.1, verbose = FALSE)
  
  # ============================================================
  # 3. 最优 pK（paramSweep）
  # ============================================================
  sweep_list  <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn       <- find.pK(sweep_stats)
  
  # 关键修复：pK 是因子，必须先转 character 再转 numeric
  optimal.pk <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(pK) %>%
    as.character() %>%
    as.numeric()
  cat("最优 pK:", optimal.pk, "\n")
  
  # ============================================================
  # 4. 同型 doublet 校正
  # ============================================================
  annotations    <- sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi       <- round(multiplet_rate * nrow(sample@meta.data))
  nExp.poi.adj   <- round(nExp.poi * (1 - homotypic.prop))
  nExp.poi.adj   <- max(nExp.poi.adj, 1)  # 至少1个
  cat("nExp:", nExp.poi, "| 校正后 nExp:", nExp.poi.adj, "\n")
  
  # ============================================================
  # 5. 计算 pANN（BiocNeighbors，修复 xtfrm bug）
  # ============================================================
  real.cells   <- rownames(sample@meta.data)
  counts       <- LayerData(sample, assay = "RNA", layer = "counts")
  data_mat     <- counts[, real.cells]
  n_real.cells <- length(real.cells)
  n_doublets   <- round(n_real.cells / (1 - 0.25) - n_real.cells)
  
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets    <- (data_mat[, real.cells1] + data_mat[, real.cells2]) / 2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data_mat, doublets)
  
  seu_wd <- CreateSeuratObject(counts = data_wdoublets)
  seu_wd <- NormalizeData(seu_wd, verbose = FALSE)
  seu_wd <- FindVariableFeatures(seu_wd, verbose = FALSE)
  seu_wd <- ScaleData(seu_wd, verbose = FALSE)
  seu_wd <- RunPCA(seu_wd, npcs = min_pc, verbose = FALSE)
  
  pca.coord <- as.matrix(seu_wd@reductions$pca@cell.embeddings[, 1:min_pc])
  nCells_wd <- nrow(pca.coord)
  rm(seu_wd); gc()
  
  # 关键修复：k 不能超过总细胞数
  k <- round(nCells_wd * optimal.pk)
  k <- max(k, 5)                          # 至少5个邻居
  k <- min(k, floor(nCells_wd / 2) - 1)  # 不超过一半细胞数
  cat("k =", k, "\n")
  
  ann_result <- BiocNeighbors::queryKNN(
    X       = pca.coord,
    query   = pca.coord[1:n_real.cells, , drop = FALSE],
    k       = k + 1,
    BNPARAM = BiocNeighbors::AnnoyParam(),
    BPPARAM = BiocParallel::MulticoreParam(workers = 4)
  )
  
  pANN_vec <- apply(
    ann_result$index[, 2:(k+1), drop = FALSE], 1,
    function(neighbors) length(which(neighbors > n_real.cells)) / k
  )
  
  # ============================================================
  # 6. 分类
  # ============================================================
  classifications <- rep("Singlet", n_real.cells)
  classifications[order(pANN_vec, decreasing = TRUE)[1:nExp.poi.adj]] <- "Doublet"
  
  cat("Doublet:", sum(classifications == "Doublet"),
      "| Singlet:", sum(classifications == "Singlet"), "\n")
  
  # ============================================================
  # 7. 返回结果
  # ============================================================
  result_df <- data.frame(
    row_names      = real.cells,
    pANN           = pANN_vec,
    doublet_finder = classifications,
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- real.cells
  return(result_df)
}

# ============================================================
# 按 sample 逐个运行
# ============================================================
# Seurat v5 修复
seu
seu <- JoinLayers(seu)

samples     <- unique(seu@meta.data$sample)
result_list <- list()

cat("In total", length(samples), "sample\n")
print(table(seu@meta.data$sample))

for (s in samples) {
  cat("\n========================================\n")
  cat("process sample:", s, "\n")
  cat("========================================\n")
  
  seu_sub <- subset(seu, subset = sample == s)
  
  tryCatch({
    result_list[[s]] <- run_doubletfinder_custom(seu_sub)
  }, error = function(e) {
    cat("ERROR in sample", s, ":", conditionMessage(e), "\n")
    result_list[[s]] <<- data.frame(
      row_names      = colnames(seu_sub),
      pANN           = NA,
      doublet_finder = "Singlet",
      stringsAsFactors = FALSE
    )
  })
}

# ============================================================
# return seu
# ============================================================
all_results           <- do.call(rbind, result_list)
rownames(all_results) <- all_results$row_names

# 按 cell name 对齐写回
seu_combined <- seu
seu_combined$pANN           <- all_results[colnames(seu), "pANN"]
seu_combined$doublet_finder <- all_results[colnames(seu), "doublet_finder"]

cat("\n===== 全部样本 Doublet 汇总 =====\n")
print(table(seu_combined$doublet_finder))
print(round(prop.table(table(seu_combined$doublet_finder)) * 100, 2))


# 查看哪些cluster doublet富集
doublet_summary <- seu_combined@meta.data %>%
  group_by(sample, doublet_finder) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = doublet_finder, values_from = n, values_fill = 0) %>%
  mutate(doublet_pct = round(Doublet / (Doublet + Singlet) * 100, 1)) %>%
  arrange(desc(doublet_pct))
doublet_summary
write.csv(doublet_summary, "doublet_summary.csv")


################################## pANN violin plot
pann_col <- grep("pANN", colnames(seu_combined@meta.data), value = TRUE)[1]

plot_df <- data.frame(
  pANN           = as.numeric(seu_combined@meta.data[[pann_col]]),
  Classification = seu_combined@meta.data$doublet_finder
)

threshold <- min(plot_df$pANN[plot_df$Classification == "Doublet"])

pdf("pANN_violin_separation.pdf", width = 3.7, height = 4)
ggplot(plot_df, aes(x = Classification, y = pANN, fill = Classification)) +
  geom_violin(alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = threshold, linetype = "dashed", 
             color = "black", linewidth = 0.8) +
  annotate("text", x = -Inf, y = threshold-0.02,
           label = round(threshold, 3),
           hjust = -0.1, vjust = 0.5)+
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  labs(title = "DoubletFinder—pANN Score Distribution",
       x = NULL, y = "pANN") +
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=12,color="black"),
    axis.text.y = element_text(size = 11, color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="top",
    legend.key.height = unit(0.6,"cm")
  )
dev.off()

################################## pANN density plot
pdf("pANN_density_separation.pdf", width = 3.7, height = 4)
ggplot(plot_df, aes(x = pANN, fill = Classification, color = Classification)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = threshold, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = threshold - 0.05, y = Inf, vjust = 1.5,
           label = paste0("threshold\n= ", round(threshold, 3)), size = 3.5) +
  scale_fill_manual(values  = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  scale_color_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  labs(title = "DoubletFinder—pANN Score Distribution",
       x = "pANN", y = "Density", fill = NULL, color = NULL) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=12,color="black"),
    axis.text.y = element_text(size = 11, color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="top",
    legend.key.height = unit(0.6,"cm")
  )
dev.off()


################################## UMAP
# 验证 doublet 列存在
print(table(seu_combined@meta.data$doublet_finder))
# 
# # ============================================================
# # 预处理 + Harmony整合 + UMAP
# # ============================================================
seu_combined <- JoinLayers(seu_combined)
seu_combined <- NormalizeData(seu_combined, verbose = FALSE)
seu_combined <- FindVariableFeatures(seu_combined, nfeatures = 2000, verbose = FALSE)
seu_combined <- ScaleData(seu_combined, verbose = FALSE)
seu_combined <- RunPCA(seu_combined, npcs = 30, verbose = FALSE)


# Harmony 整合批次效应
table(seu_combined@meta.data$sample)
seu_combined <- RunHarmony(seu_combined, group.by.vars = "sample", verbose = FALSE)
DefaultDimReduc(seu_combined)
# UMAP
seu_combined
seu_combined <- RunUMAP(seu_combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

# ============================================================
# ============================================================
seu_combined
umap_df <- as.data.frame(seu_combined@reductions$umap@cell.embeddings)

colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$doublet_finder <- seu_combined@meta.data$doublet_finder
umap_df$sample         <- seu_combined@meta.data$sample

df_singlet <- dplyr::filter(umap_df, doublet_finder == "Singlet")
df_doublet <- dplyr::filter(umap_df, doublet_finder == "Doublet")

# 图1：总览
pdf("UMAP_doublet.pdf", width = 3.4, height = 4)
ggplot() +
  geom_point(data = df_singlet, aes(x = UMAP_1, y = UMAP_2, color = doublet_finder),
             size = 0.1, alpha = 0.4) +
  geom_point(data = df_doublet, aes(x = UMAP_1, y = UMAP_2, color = doublet_finder),
             size = 0.0000001, alpha = 0.9) +
  scale_color_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c"),
                     labels = c( paste0("Doublet (n=", nrow(df_doublet), ")"),
                                 paste0("Singlet (n=", nrow(df_singlet), ")"))) +
  labs(title = "DoubletFinder — UMAP", x = "UMAP 1", y = "UMAP 2", color = NULL) +
  theme_bw() +
  theme(
    panel.border    = element_rect(fill = NA, color = "black", linewidth = 1.2),
    plot.title      = element_text(size = 13, hjust = 0.5),
    axis.text       = element_text(size = 10, color = "black"),
    panel.grid      = element_blank(),
    legend.position = "top",
    legend.key      = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
dev.off()


################################## UMAP split
# 2. by sample
pdf("UMAP_doublet_split.pdf", width = 11.5, height = 4)
ggplot() +
  geom_point(data = df_singlet, aes(x = UMAP_1, y = UMAP_2, color = doublet_finder),
             size = 0.1, alpha = 0.3) +
  geom_point(data = df_doublet, aes(x = UMAP_1, y = UMAP_2, color = doublet_finder),
             size = 0.8, alpha = 0.9) +
  scale_color_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  facet_wrap(~ sample, ncol = 4) +
  labs(title = "DoubletFinder — UMAP by sample",
       x = "UMAP 1", y = "UMAP 2", color = NULL) +
  theme_bw() +
  theme(
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 1.2),
    plot.title       = element_text(size = 13, hjust = 0.5),
    axis.text        = element_text(size = 8, color = "black"),
    panel.grid       = element_blank(),
    legend.position  = "top",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text       = element_text(size = 9, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
dev.off()





################################## quality control comparison
# Check how doublets singlets differ in QC measures per sample.
seu_combined <- PercentageFeatureSet(seu_combined, pattern = "^MT.", col.name = "percent.mt")
pdf('QC_doubletVSsinglet.pdf', width = 12, height = 4)
VlnPlot(
  seu_combined, 
  group.by = 'sample', 
  split.by = "doublet_finder",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3, 
  pt.size = 0,
  cols = c("#e31a1c", "#1f78b4")
) +
  theme(legend.position = "right",
        plot.title = element_text(size=13,hjust = 0.5),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(size=12, angle = 70, hjust = 1,color="black"),
        axis.text.y = element_text(size = 11, color="black"),
        legend.title = element_blank())
dev.off()




################################## doublet number summary
# Get doublets per sample
doublets_summary <- seu_combined@meta.data %>% 
  group_by(sample, doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(sample) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
write.csv(doublets_summary, "doublets_summary-before.csv")
################################### stack  ###################################
data <- read.csv("doublets_summary-before.csv")[,-1]
head(data)


pdf('doublet_finder_percent.pdf', width = 3.1, height = 3)
data$percent <- as.numeric(sub("%","", data$percent))
ggplot(data, aes(x = sample, y = percent, fill = doublet_finder)) +
  geom_bar(stat = "identity" , alpha = 0.7) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  ylab("Percent of Cells") +
  xlab("") +
  # ggtitle("Stacked Doublet/Singlet Count per sample") +
  ggtitle("")+
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=12, angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(size = 11, color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="right",
    legend.key.height = unit(0.6,"cm")
  )
dev.off()


pdf('doublet_finder_number.pdf', width = 2.1, height = 3.5)
ggplot(data, aes(x = sample, y = total_count, fill = doublet_finder)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  ylab("Number of Cells") +
  xlab("") +
  # ggtitle("Stacked Doublet/Singlet Count per sample") +
  ggtitle("")+
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=12, angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(size = 11, color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="top",
    legend.key.height = unit(0.6,"cm")
  )
dev.off()









############################################# violin plot  #############################################
meta <- seu_combined@meta.data
head(meta)

summary_df_filter <- meta%>%
  group_by(sample,doublet_finder) %>%
  summarise(
    nCount_q5 = quantile(nCount_RNA, 0.05, na.rm = TRUE),
    nCount_q95 = quantile(nCount_RNA, 0.95, na.rm = TRUE),
    
    nFeature_q5 = quantile(nFeature_RNA, 0.05, na.rm = TRUE),
    nFeature_q95 = quantile(nFeature_RNA, 0.95, na.rm = TRUE),
    
    percent_mt_q5 = quantile(percent.mt, 0.05, na.rm = TRUE),
    percent_mt_q95 = quantile(percent.mt, 0.95, na.rm = TRUE)
  )
write.csv(summary_df_filter,"summary_median_doublet.csv")




pdf("nCount_RNA_doubletfinder.pdf", width = 4, height = 3.2)
p1<- ggplot(meta, aes(x = sample, y = nCount_RNA, fill = doublet_finder)) +
  geom_violin(
    position = position_dodge(width = 0.8),
    scale = "width",
    trim = TRUE,
    width = 0.7,
    alpha = 0.7, 
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    position = position_dodge(width = 0.8),
    size = 0.2
  )+
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
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


head(meta)
pdf("nFeature_RNA_doubletfinder.pdf", width = 4, height = 3.2)
p2<- ggplot(meta, aes(x = sample, y = nFeature_RNA, fill = doublet_finder)) +
  geom_violin(
    position = position_dodge(width = 0.8),
    scale = "width",
    trim = TRUE,
    width = 0.7,
    alpha = 0.7, 
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    position = position_dodge(width = 0.8),
    size = 0.2
  )+
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
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





head(meta)
pdf("percent.mt_doubletfinder.pdf", width = 4, height = 3.2)
p3 <- ggplot(meta, aes(x = sample, y = percent.mt, fill = doublet_finder)) +
  geom_violin(
    position = position_dodge(width = 0.8),
    scale = "width",
    trim = TRUE,
    width = 0.7,
    alpha = 0.7, 
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    position = position_dodge(width = 0.8),
    size = 0.2
  )+
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
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



pdf("QC_doubletfinder.pdf", width = 10, height = 3)
# 组合图
combined <- ggarrange(
  p1, p2, p3,
  widths = c(1, 1, 0.96),
  ncol = 3, nrow = 1
)
combined
dev.off()








############################################# final  #############################################
data <- read.csv("doublets_summary-before.csv")[,-1]
head(data)


p <- ggplot(data, aes(x = sample, y = total_count, fill = doublet_finder)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("Singlet" = "#1f78b4", "Doublet" = "#e31a1c")) +
  ylab("Number of Cells") +
  xlab("") +
  # ggtitle("Stacked Doublet/Singlet Count per sample") +
  ggtitle("")+
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=12, angle = 45, hjust = 1,color="black"),
    axis.text.y = element_text(size = 11, color="black"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position = "none",   # x = 0.5 (水平居中), y = 1.2 (往上移)
    legend.justification = c(0.5, 0),
    legend.key.height = unit(0.6,"cm"),
    # plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
  )
p



pdf("QC_doubletfinder_final.pdf", width = 12, height = 3.5)
# 组合图
combined <- ggarrange(
  p, p1, p2, p3,
  widths = c(0.75, 1, 1, 0.95),
  ncol = 4, nrow = 1
)
combined
dev.off()


# ============================================================
# remove doublet and save singlet
# ============================================================
meta <- seu_combined@meta.data
write.csv(meta,"meta_doublefinder.csv")

save(seu_combined,file = "/Users/birzh586/Desktop/skin-horse/data/seu_with_doublet_labels.rds")

seu_singlet <- subset(seu_combined, subset = doublet_finder == "Singlet")
cat("post filter：", ncol(seu_singlet), "\n")

seu_singlet
# 提取 counts
counts <- GetAssayData(seu_singlet, layer = "counts")

# 重新创建 Seurat object
seu_singlet_clean <- CreateSeuratObject(
  counts = counts,
  meta.data = seu_singlet@meta.data,
  project = "seu_singlet"
)

# 查看结果
seu_singlet_clean
# An object of class Seurat 
# 30152 features across 85574 samples within 1 assay 
# Active assay: RNA (30152 features, 0 variable features)
# 1 layer present: counts
saveRDS(seu_singlet_clean, file = "/Users/birzh586/Desktop/skin-horse/data/seu_singlet_clean.rds")



