library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(Seurat)


setwd("/Users/birzh586/Desktop/skin-horse/analysis")
dir.create("3.Cell-annotation")
setwd("/Users/birzh586/Desktop/skin-horse/analysis/3.Cell-annotation")
merge_Singlet <- readRDS("/Users/birzh586/Desktop/skin-horse/data/seu_singlet_clean.rds")
merge_Singlet
# An object of class Seurat 
# 30152 features across 85574 samples within 1 assay 
# Active assay: RNA (30152 features, 0 variable features)
# 1 layer present: counts

# 去除核糖体和线粒体基因
genes_to_remove <- c(grep("^RPS", rownames(merge_Singlet), value = TRUE),
                     grep("^RPL", rownames(merge_Singlet), value = TRUE),
                     grep("^MT-", rownames(merge_Singlet), value = TRUE))

merge <- merge_Singlet[!rownames(merge_Singlet) %in% genes_to_remove, ]
merge
# An object of class Seurat 
# 30058 features across 85574 samples within 1 assay 
# Active assay: RNA (30058 features, 0 variable features)
# 1 layer present: counts

############################### merge ###############################
## normalise and scale
merge <- NormalizeData(merge)

## Calculate cell-cycle scores
merge <- CellCycleScoring(object = merge, 
                          g2m.features = cc.genes$g2m.genes,
                          s.features = cc.genes$s.genes)

pdf("CellCycleScoring.pdf", width = 6.6, height = 4)
VlnPlot(merge, features = c("S.Score", "G2M.Score"), group.by = "sample", ncol = 2, pt.size = 0.1, raster = T)
dev.off()

table(merge$Phase, merge$sample)
#      C_1a  C_1b  C_2a  C_2b
# G1   5379 16251 10661 17198
# G2M  1970  4436  4892  3911
# S    2253  6507  6073  6043
## Feature selection
## merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 3000, verbose = FALSE, assay = "RNA")
## top20 <- head(VariableFeatures(merge), 20)
## LabelPoints(plot = VariableFeaturePlot(merge), points = top20, repel = TRUE)


## Step 1: Automatic selection of high variable genes
merge <- FindVariableFeatures(
  merge,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE,
  assay = "RNA"
)
current_features <- VariableFeatures(merge)


## Step 2: Manually add genes of interest
manual_genes <- c("FCER1A", "MS4A2", "HPGDS", "LTC4S", "IL1RL1", # mast cell
                 "CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD7")  # T cell  

# Step 3: Merge and remove duplicates
new_features <- unique(c(current_features, manual_genes))

# update VariableFeatures
VariableFeatures(merge) <- new_features

## Step 3: 确认添加成功
manual_genes %in% VariableFeatures(merge)


# Scale: Z-score transformation
merge <- ScaleData(merge, vars.to.regress = c("percent.mt", "nCount_RNA"),assay = "RNA")

## RunPCA
merge <- RunPCA(merge, features = new_features, assay = "RNA", verbose = FALSE)

## RunHarmony
# merge <- RunHarmony(merge, group.by.vars = "sample", max.iter.harmony = 10)
head(merge@meta.data)
merge <- RunHarmony(merge, 
                    group.by.vars = "sample", 
                    max.iter.harmony = 10,  # Increase iterations
                    theta = 2,              # Increase diversity penalty
                    lambda = 1,             # Adjust ridge regression penalty
                    sigma = 0.1,            # Adjust soft clustering
                    verbose = TRUE)


ElbowPlot(merge, reduction = "harmony", ndims = 50)
cumsum(merge@reductions$harmony@stdev)  

sum(merge@reductions$harmony@stdev[1:10])/sum(merge@reductions$harmony@stdev)
sum(merge@reductions$harmony@stdev[1:20])/sum(merge@reductions$harmony@stdev)
sum(merge@reductions$harmony@stdev[1:30])/sum(merge@reductions$harmony@stdev)
sum(merge@reductions$harmony@stdev[1:40])/sum(merge@reductions$harmony@stdev)
sum(merge@reductions$harmony@stdev[1:50])/sum(merge@reductions$harmony@stdev)

# RunUMAP
merge <- RunUMAP(merge, 
                 reduction = 'harmony', 
                 dims = 1:30,
                 n.neighbors = 30,      # Increase for smoother manifold
                 min.dist = 0.3,        # Adjust for cluster spacing
                 metric = "cosine",
                 verbose = TRUE)


set.seed(1)
merge <- FindNeighbors(merge, 
                       reduction = "harmony", 
                       dims = 1:30,
                       k.param = 20)       # Adjust k-nearest neighbors

set.seed(1)
merge <- FindClusters(merge, 
                      resolution = 0.8,    # Try higher resolution first
                      algorithm = 1)       # Louvain algorithm

table(merge$seurat_clusters)



DimPlot(merge, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
FeaturePlot(merge, features = c("CD3D", "CD3E", "MS4A2", "FCER1A", "AIF1", "CTSS"), reduction = "umap")

# Umap
pdf("umap_cluster.pdf", width = 6, height = 5)
DimPlot(merge, reduction = "umap",pt.size = 1,  label = TRUE, raster= T)
dev.off()



# Cell Type Assignation
library(fgsea)
hsaPanglaoDB <- gmtPathways('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/hsaPanglaoDB.gmt')

cellTypeAssignation <- lapply((seq_along(unique(Idents(merge))) - 1), function(cID){
  selectedCells <- Idents(merge) %in% cID
  otherAvg <- rowMeans(merge@assays$RNA$counts[,!selectedCells] != 0)
  cellAvg <- rowMeans(merge@assays$RNA$counts[,selectedCells] != 0)
  log2Freq <- log2(cellAvg/otherAvg)
  log2Freq <- log2Freq[is.finite(log2Freq)]
  E <- fgseaMultilevel(hsaPanglaoDB, log2Freq, eps = 0)
  E <- E[E$NES > 0,]
  E <- E[order(E$padj),]
})

new.cluster.ids <- unlist(lapply(cellTypeAssignation, function(x) {
  as.data.frame(x)[1, 1]
}))


## cluster.ids <- c("Hepatocytes", new.cluster.ids[-1])
head(cellTypeAssignation)

cellTypeAssignation[[1]][1] # Endothelial cells
cellTypeAssignation[[2]][1] # Enteroendocrine cells
cellTypeAssignation[[3]][1] # Keratinocytes
cellTypeAssignation[[4]][1] # Noradrenergic neurons
cellTypeAssignation[[5]][1] # Myoepithelial cells
cellTypeAssignation[[6]][1] # Keratinocytes
cellTypeAssignation[[7]][1] # 
cellTypeAssignation[[8]][1] # 
cellTypeAssignation[[9]][1] # 
cellTypeAssignation[[10]][1] # 
cellTypeAssignation[[11]][1] # 
cellTypeAssignation[[12]][1] # Smooth muscle cells
cellTypeAssignation[[13]][1] # 
cellTypeAssignation[[14]][1] # Fibroblasts
cellTypeAssignation[[15]][1] # 
cellTypeAssignation[[16]][1] # Melanocytes
cellTypeAssignation[[17]][1] # Mast cells
cellTypeAssignation[[18]][1] # T cells
cellTypeAssignation[[19]][1] # Eosinophils
cellTypeAssignation[[20]][1] # 
cellTypeAssignation[[21]][1] # Macrophages
cellTypeAssignation[[22]][1] # 

############### Cell marker genes
# Compute differentiall expression
merge <- JoinLayers(merge)
DefaultAssay(merge)        # 是否是 RNA
# Idents(merge)              # cluster 是否存在
Layers(merge)              # 是否 >1
markers_genes <- FindAllMarkers(merge, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")

write.csv(markers_genes,"markers_genes.csv")

# We can now select the top 25 up regulated genes for plotting.
markers_genes %>%
  group_by(cluster) %>%
  top_n(-25, p_val_adj) -> top25
top25

# now select the top 25 up regulated genes for plotting
# mypar(2, 5, mar = c(4, 6, 3, 1))
# for (i in unique(top25$cluster)) {
#   barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
#           horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
#   abline(v = c(0, 0.25), lty = c(1, 2))
# }


# visualize them as a heatmap
markers_genes %>%
  group_by(cluster) %>%
  top_n(-5, p_val_adj) -> top5

# create a scale.data slot for the selected genes
# merge <- ScaleData(merge, features = as.character(unique(top5$gene)), assay = "RNA")
# pdf('DEG_heatmap.pdf', width = 23, height = 12)
# DoHeatmap(merge, features = as.character(unique(top5$gene)), assay = "RNA")
# dev.off()




####################################### cell annotation  ####################################### 
# Basal Keratinocytes
pdf('Basal_Keratinocytes_markers.pdf', width = 16, height = 9)
basal_markers <- c("KRT5","KRT14","KRT15","COL17A1","ITGB4","ITGA6","TP63","LAMB3","LAMA3","DST","PLEC","POSTN","KLF5","CD44")
FeaturePlot(merge,  features = basal_markers, raster = T) 
dev.off()


# Stress_response Keratinocytes
pdf('Stress_response_markers.pdf', width = 16, height = 9)
Stress_response_markers <- c("FOS", "JUN", "JUNB", "ATF3", "EGR1", "DUSP1", "MYC", "NFKBIA", "GADD45B", "KLF4",
                       "CCL2", "CCL27", "IL1R2","HSPA1B" )
FeaturePlot(merge,  features = Stress_response_markers, raster = T) 
dev.off()


# Proliferative Keratinocytes
pdf('Proliferative Keratinocytes.pdf', width = 16, height = 9)
proliferative_markers <- c("MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1", "CDK1", "TYMS", "RRM2", "CENPA", "NUSAP1")
FeaturePlot(merge,  features = proliferative_markers, raster = T) 
dev.off()

# hyperproliferative Keratinocytes
pdf('Hyperproliferative Keratinocytes.pdf', width = 16, height = 9)
hyperproliferative_markers <- c("KRT6A","KRT6B","KRT6C", "KRT16", "KRT17")
FeaturePlot(merge,  features = hyperproliferative_markers, raster = T) 
dev.off()




# Endothelial Cells  "CCN1", "CCN2", "SOX18","EGFL7", "SOX18", "S1PR1"
pdf('Blood_Endothelial_markers.pdf', width = 16, height = 9)
endothelial_markers <- c("PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "PLVAP", "CD34", "EMCN", "ERG", "ESAM")
FeaturePlot(merge,  features = endothelial_markers, raster = T) 
dev.off()

# lymphatic endothelial cells (LECs).
pdf('Lymphatic_endothelial_markers.pdf', width = 16, height = 9)
Lymphatic_endothelial_markers <- c("PROX1", "LYVE1", "FLT4",  "STAB2")
FeaturePlot(merge,  features = Lymphatic_endothelial_markers, raster = T) 
dev.off()




# Basal Progenitor Keratinocytes
pdf('Basal Progenitor Keratinocytes.pdf', width = 8, height = 6.3)
Basal_Progenitor_markers  <- c("COL17A1","ITGB4","LGR5","THY1", "BCAM", "BGN", "POSTN","SPARC")
FeaturePlot(merge,  features = Basal_Progenitor_markers, raster = T) 
dev.off()
# Basal Progenitor-like Keratinocytes
pdf('Basal Progenitor-like  Keratinocytes.pdf', width = 8, height = 6.3)
Basal_Progenitor_like_markers <- c("KRT8","WNT6","LBH","PPARG","IL1R2","DSG2")
FeaturePlot(merge,  features = Basal_Progenitor_like_markers, raster = T) 
dev.off()
# Basement Membrane-associated Keratinocytes/ECM-interacting keratinocytes
pdf('Basement Membrane-associated Keratinocytes.pdf', width = 12, height = 6.3)
Basement_Membrane_markers <- c("COL17A1", "LAMB3", "LAMA3", "ITGB4", "KRT15", "SOX9","POSTN", "SPARC", "BGN", "ECM1") 
FeaturePlot(merge,  features = Basement_Membrane_markers, raster = T) 
dev.off()
# ECM_remodeling keratinocytes 
pdf('ECM_remodeling keratinocytes.pdf', width = 8, height = 6.3)
ECM_remodeling_markers <- c("POSTN", "SPARC", "SPARCL1", "CCN1","CCN2", "BGN", "SOX9", "MYL9")
FeaturePlot(merge,  features = ECM_remodeling_markers, raster = T) 
dev.off()




# Sebocytes / Sebaceous Gland Cells 皮脂腺细胞 /Sebaceous Gland Sebocytes (皮脂腺皮脂细胞)
pdf('Sebocytes.pdf', width = 16, height = 9)
sebocyte_markers <- c("SCD", "THRSP", "ELOVL3", "AWAT1", "AWAT2", "SOAT1", "CIDEA", "PGRMC1", "FASN", "PPARG","KRT7")
FeaturePlot(merge,  features = sebocyte_markers, raster = T) 
dev.off()
# Sweat Gland Cells 汗腺细胞
pdf('Sweat Gland Cells.pdf', width = 16, height = 9)
Sweat_Gland_markers <- c("PIP", "AZGP1", "CRISP3","AQP5", "ATP1A1", "ATP1B1", "SLC12A2", "EHF","KRT7", "KRT8", "KRT18", "KRT19", "AQP5", "CA8", "AZGP1", "CLDN3" )
FeaturePlot(merge,  features = Sweat_Gland_markers, raster = T) 
dev.off()
# Eccrine Sweat Gland Cells 小汗腺细胞
pdf('Eccrine Sweat Gland Cells.pdf', width = 16, height = 9)
eccrine_markers <- c("AQP5", "ATP1A1", "ATP1B1", "SLC12A2", "EHF") # "CLDN8", "ELF3", "KCNE2","SCNN1B"
FeaturePlot(merge,  features = eccrine_markers, raster = T) 
dev.off()
# Apocrine Sweat Gland Epithelial Cells 大汗腺上皮细胞
pdf('Apocrine Sweat Gland Cells.pdf', width = 16, height = 9)
apocrine_markers <- c("PIP", "AZGP1", "CRISP3" )
FeaturePlot(merge,  features = apocrine_markers, raster = T) 
dev.off()

gland_comparison <- c(
  # 小汗腺特异
  "EHF", "CLDN8", "ATP1A1",
  # 大汗腺特异  
  "PIP", "SCGB2A2", "AZGP1",
  # 皮脂腺特异
  "SCD", "AWAT1", "CIDEA"
)
FeaturePlot(merge, features = gland_comparison, raster = T)






# Hair shaft keratinocytes (cortex/cuticle lineage)
pdf('Hair shaft keratinocytes.pdf', width = 8, height = 6.3)
Hair_shaft_markers <- c("KRT31", "KRT32", "KRT35", "KRT85", "KRT89", "HOXC13",  "DSG4") # "DLX3",
FeaturePlot(merge,  features = Hair_shaft_markers, raster = T) 
dev.off()
# Outer Root Sheath
pdf('Outer Root Sheath.pdf', width = 16, height = 9)
ors_markers <- c("KRT71", "KRT72", "KRT73", "KRT25", "GATA3", "TCHH", "IVL", "DLX3", "GSDMA", "KRT27")
FeaturePlot(merge,  features = ors_markers, raster = T) 
dev.off()
# Inner Root Sheath
pdf('Inner Root Sheath.pdf', width = 16, height = 9)
irs_markers <- c("KRT23", "KRT75",  "FOXQ1", "GRHL1", "TCHH", "DSG4", "PADI3", "SERPINB13", "IVL")
FeaturePlot(merge,  features = irs_markers, raster = T) 
dev.off()






# Terminally differentiated / cornified keratinocytes
pdf('late_differentiated_KCs.pdf', width = 10, height = 6.3)
# late_diff_KC_markers <- c("S100A7","S100A12", "KLK8",  "CRABP2",  "BBOX1","KLK13","FLG","LOR")

late_diff_KC_markers <-c(
  "FLG",
  "LOR",
  "IVL",
  "SPRR1A", "SPRR1B",
  "SPRR2A", "SPRR2B", "SPRR2C", "SPRR2D", "SPRR2E", "SPRR2F", "SPRR2G",
  "SPRR3",
  "TGM1", "TGM3",
  "KRT9",
  "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F",
  "LCE2A", "LCE2B", "LCE2C", "LCE2D",
  "LCE3A", "LCE3B", "LCE3C", "LCE3D", "LCE3E",
  "LCE4A", "LCE5A"
)
FeaturePlot(merge,  features = late_diff_KC_markers, raster = T) 
dev.off()



# intermediate_differentiated_KCs
pdf('intermediate_differentiated_KCs.pdf', width = 10, height = 6.3)
intermediate_diff_KC_markers <- c("SBSN", "SFN", "SPINK5",  "KLK7",
  "CASP14",
  "AQP3",
  "ELOVL1",
  "ELOVL3",
  "HACD3",
  "PRDX5",
  "GPX3",
  "LYPD3",
  "TACSTD2",
  "MSMB",
  "DMKN",
  "EMP2",
  "GLTP",
  "PDZK1IP1"
)

FeaturePlot(merge,  features = intermediate_diff_KC_markers, raster = T) 
dev.off()


pdf('early_differentiated_KCs.pdf', width = 10, height = 6.3)
early_diff_KC_markers <- c("KRT1","KRT2","DSC1","DSG1","KLK5","CASP14","CALML5","LY6D")
FeaturePlot(merge,  features = early_diff_KC_markers, raster = T) 
dev.off()





# inflammatory_spinous_KCs
pdf('inflammatory_spinous_KCs.pdf', width = 7, height = 8.8)
inflammatory_spinous_KCs <- c("S100A7","S100A12","CRABP2","BPIFC","PGLYRP3","FOXQ1")
FeaturePlot(merge,  features = inflammatory_spinous_KCs, raster = T) 
dev.off()
# activated_spinous_KCs
pdf('activated_spinous_KCs.pdf', width = 10, height = 6.3)
activated_spinous_KCs <- c("ATP6V1C2","UNC93A","KRT6B","KRT17","PLET1","TNFAIP8L3","GRHL1")
FeaturePlot(merge,  features = activated_spinous_KCs, raster = T) 
dev.off()
# spinous_KCs
pdf('spinous_KCs.pdf', width = 7, height = 8.8)
spinous_KCs <- c("KRT1","KRT2B","KRT77","DSC1","CALML5","SBSN")
FeaturePlot(merge,  features = spinous_KCs, raster = T) 
dev.off()






# suprabasal_markers  (基底上层角质形成细胞) 所有离开基底层的角质形成细胞，包括: 棘层 + 颗粒层 + 角质层
pdf('suprabasal KCs.pdf', width = 8, height = 6.3)
suprabasal_markers <- c(
  # 分化角蛋白（离开基底层的标志）
  "KRT1", "KRT10", "KRT2",
  # 早期分化marker
  "IVL", "DSG1", "DSC1",
  # 代谢转换
  "TGM1", "ALOX12B",
  # vs 基底层（应该低表达或无）
  "KRT5", "KRT14", "TP63" 
)
FeaturePlot(merge,  features = suprabasal_markers, raster = T) 
dev.off()

# Granular Keratinocytes 颗粒层角质形成细胞 位置: 棘层之上，角质层之下 角化颗粒(keratohyalin granules) - FLG, LOR
pdf('Granular Keratinocytes.pdf', width = 13.5, height = 9)
granular_markers <- c(
  # 角化颗粒核心蛋白（金标准）
  "FLG",      # Filaggrin - 最重要
  "LOR",      # Loricrin - 最丰富
  "FLG2",     # Filaggrin-2
  # 其他颗粒层标志
  "KRTDAP",   # 角蛋白相关蛋白
  "HRNR",     # Hornerin
  # 蛋白酶
  "CASP14",   # 切割profilaggrin
  "KLK7",     # 角质层脱落
  # 脂质屏障
  "ABCA12",   # 脂质转运
  "ALOX12B",  # 脂氧合酶
  # 其他
  "CSTA", "SPINK5", "SBSN"
)
FeaturePlot(merge,  features = granular_markers, raster = T) 
dev.off()

# Spinous Keratinocytes  (棘层角质形成细胞) 位置: 基底层之上，颗粒层之下
pdf('Spinous Keratinocytes.pdf', width = 16, height = 9)
spinous_markers <- c(
  # 棘层角蛋白
  "KRT1", "KRT10",
  # 桥粒蛋白（棘层丰富）
  "DSG3", "DSC3",  # 下棘层
  "DSG1", "DSC2",  # 上棘层
  "DSP", "PKP1", "JUP",
  # 缝隙连接
  "GJB2", "GJB6",
  # 早期分化marker
  "IVL",  # 在上棘层开始表达
  "SPRR1A", "SPRR1B"
)

spinous_markers <- c("PLET1","SPINK5","CLDN1","SBSN","GRHL1","SCEL","LAD1")
FeaturePlot(merge,  features = spinous_markers, raster = T) 
dev.off()

# Differentiated Keratinocytes
pdf('Differentiated_markers.pdf', width = 16, height = 9)
differentiated_markers <- c(
  # 分化角蛋白（核心标志）
  "KRT1", "KRT10",
  # 晚期分化角蛋白
  "KRT2", "KRT9",
  # 角化包膜蛋白（分化进程）
  "IVL",      # 早期
  "LOR",      # 晚期
  "SPRR3",    # 中期
  # 桥粒转换（基底→分化）
  "DSG1", "DSC1",  # 分化型
  # 脂质代谢（屏障形成）
  "ALOX12B", "ALOXE3",
  # 蛋白酶
  "KLK5", "KLK7", "CASP14"
)
FeaturePlot(merge,  features = differentiated_markers, raster = T) 
dev.off()

#  DotPlot
pdf('Keratinocyte_Differentiation_DotPlot.pdf', width = 14, height = 8)
all_kc_markers <- c(
  # 基底层
  "KRT5", "KRT14", "TP63",
  # 棘层
  "KRT1", "KRT10", "DSG3", "IVL",
  # 颗粒层
  "FLG", "LOR", "CASP14",
  # 晚期
  "SBSN", "KRT2", "SCEL"
)
all_kc_markers_filtered <- all_kc_markers[all_kc_markers %in% rownames(merge)]
DotPlot(merge, features = all_kc_markers_filtered, 
        cols = c("lightgrey", "red"),
        dot.scale = 8) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ===== 5. 分化梯度对比 =====
pdf('Differentiation_Gradient.pdf', width = 20, height = 12)
gradient_markers <- c(
  # 基底层
  "KRT5", "KRT14", "TP63",
  # 早期分化
  "KRT1", "KRT10", "IVL",
  # 颗粒层
  "FLG", "LOR", "CASP14",
  # 晚期
  "SBSN", "SCEL", "KRT2"
)
gradient_markers_filtered <- gradient_markers[gradient_markers %in% rownames(merge)]
FeaturePlot(merge, features = gradient_markers_filtered, raster = T, ncol = 4)
dev.off()




# Smooth Muscle Cells
pdf('Smooth Muscle Cells.pdf', width = 16, height = 9)
smooth_muscle_markers <- c("ACTA2", "MYH11", "TAGLN", "TPM1", "TPM2", "CNN1", "MYLK", "MYL9", "ACTG2", "SMTN")
FeaturePlot(merge,  features = smooth_muscle_markers, raster = T) 
dev.off()
# Fibroblasts
pdf('Fibroblasts.pdf', width = 16, height = 9)
fibroblast_markers <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "DCN", "LUM", "PDGFRA", "THY1", "FAP")
FeaturePlot(merge,  features = fibroblast_markers, raster = T) 
dev.off()
# Melanocytes
pdf('Melanocytes.pdf', width = 16, height = 9)
melanocyte_markers <- c("PMEL", "MLANA", "TYRP1", "TYR", "DCT", "MITF", "SOX10", "OCA2", "RAB27A", "GPNMB")
FeaturePlot(merge,  features = melanocyte_markers, raster = T) 
dev.off()
# T Cells
pdf('T Cells.pdf', width = 16, height = 9)
t_cell_markers <- c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD7", "TRAC", "TRBC1", "LCK", "IL7R")
FeaturePlot(merge,  features = t_cell_markers, raster = T) 
dev.off()
# Macrophages
pdf('Macrophages.pdf', width = 16, height = 9)
macrophage_markers <- c("CD68", "CD163", "AIF1", "FCGR3A", "MARCO", "MSR1", "CD14", "LYZ", "TYROBP", "CTSS")
FeaturePlot(merge,  features = macrophage_markers, raster = T) 
dev.off()
# Mast cells
pdf('Mast cells.pdf', width = 8, height = 6.3)
mast_markers <- c("FCER1A", "MS4A2", "HPGDS", "LTC4S")
FeaturePlot(merge,  features = mast_markers, raster = T) 
dev.off()




###### annotation 
saveRDS(merge, file = "/Users/birzh586/Desktop/skin-horse/data/merge-umap.rds")
merge <- readRDS("/Users/birzh586/Desktop/skin-horse/data/merge-umap.rds")
merge
# An object of class Seurat 
# 30058 features across 87218 samples within 1 assay 
# Active assay: RNA (30058 features, 2003 variable features)
# 3 layers present: counts, data, scale.data
# 4 dimensional reductions calculated: pca, harmony, tsne, umap




# ## subset
# merge <- subset(x = merge, idents =c(0:7))
# # 29093 features across 57210 samples within 1 assay 
# table(merge$seurat_clusters)
# Idents(merge) <- "seurat_clusters"
# # levels(Idents(seurat.integrated))[levels(Idents(seurat.integrated)) %in% "13"] <- "12"

# merge <- subset(x = merge, idents =c(0:25))
table(Idents(merge))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22 
# 7984 7444 7098 6367 5507 5007 4997 4916 4623 4459 3808 2646 2404 2199 2191 2039 1474 1400 1369 1280 1277 1252  960 
# 23   24   25   26   27   28 
# 905  692  454  318  275  229 


merge <- RenameIdents(object=merge,
                      "2" = "Stress response KCs",
                      "20" = "Stress response KCs",
                      
                      "0" = "Hyperproliferative KCs",
                      "18" = "Hyperproliferative KCs",
                      "28" = "Hyperproliferative KCs",
                      
                      "1" = "Basal KCs",
                      "8" = "Endothelial", 
                      "10" = "Endothelial", 
                      "12" = "Endothelial", 
                      
                      "3" = "Spinous KCs", 
                      "16" = "Activated spinous KCs", 
                      "23" = "Inflammatory spinous KCs",
                      
                      "5" = "Proliferative KCs",
                      "4" = "Transitional KCs",
                      "7" = "Apocrine sweat gland cells",
                      "24" = "Eccrine sweat gland cells",
                      
                      "6" = "Hair shaft KCs",
                      "17" = "Hair shaft KCs",
                      
                      "9" = "Sebocytes",
                      "11" = "Basal progenitor-like KCs",
                      
                      "14" = "Outer root sheath KCs",
                      "26" = "Inner root sheath KCs",
                      
                      "13" = "Smooth muscle cells",
                      "27" = "Smooth muscle cells",
                      "15" = "Fibroblasts",
                      "19" = "Melanocytes",
                      "21" = "Mast cells",
                      "22" = "T cells",
                      "25" = "Macrophages"
)




# Store the current identities in a new column of meta.data called CellType
# 统计每个cluster细胞数
cluster_counts <- table(Idents(merge))

# 按数量从大到小排序
new_levels <- names(sort(cluster_counts, decreasing = TRUE))

# 重新设定顺序
merge <- SetIdent(merge, value = factor(Idents(merge), levels = new_levels))
table(merge@active.ident)
# Endothelial     Hyperproliferative KCs        Stress response KCs                  Basal KCs 
# 10835                       9582                       8375                       7444 
# Hair shaft KCs                Spinous KCs           Transitional KCs          Proliferative KCs 
# 6397                       6367                       5507                       5007 
# Apocrine sweat gland cells                  Sebocytes  Basal progenitor-like KCs        Smooth muscle cells 
# 4916                       4459                       2646                       2474 
# Outer root sheath KCs                Fibroblasts      Activated spinous KCs                Melanocytes 
# 2191                       2039                       1474                       1280 
# Mast cells                    T cells   Inflammatory spinous KCs  Eccrine sweat gland cells 
# 1252                        960                        905                        692 
# Macrophages      Inner root sheath KCs 
# 454                        318 


table <- table(merge@active.ident,merge$sample) %>% data.frame()
colnames(table) <- c("Cells", "sample","Freq")
write.csv(table, "table_sample_cell.csv")
saveRDS(merge, file = "/Users/birzh586/Desktop/skin-horse/data/merge_annotation.rds")




load("merge_annotation.RData")
######################### DimPlot
# allcolour <- c("#C06C84","#b2db87","#f7905a","#e187cb","wheat1","#a48cbe","#e2b159")
#
# allcolour <- c(
#   "#C06C84", # red
#   "#4DBBD5", # teal
#   "springgreen4",
#   "#a48cbe", # blue
#   "#F39B7F", # orange-pink
#   "#DC0000", # dark red
#   "#7E6148", # brown
#   "#e187cb", # navy
#   "#DF8F44", # amber
#   "#B24745", # brick red
#   "#91D1C2", # light teal
#   "#8491B4", # grey-blue
#   "#79AF97", # muted green
#   "#6A6599", # purple
#   "#FFDC91", # pale yellow
#   "slategrey" , # grey
#   "#1E90FF"  # bright blue
# )
#
#
# allcolour <- c(
#   "#C06C84", # 柔和红
#   "#4DBBD5", # 柔和蓝绿
#   "#4CAF50", # 绿色
#   "#a48cbe", # 蓝紫
#   "#F39B7F", # 橙粉
#   "#7E6148", # 棕色
#   "#FFD700", # 金黄
#   "darkgoldenrod", # 粉紫
#   "#DF8F44", # 琥珀
#   "#B24745", # 砖红
#   "#e187cb", # 灰蓝
#   "#79AF97", # 暗绿色
#   "#91D1C2", # 浅蓝绿
#   "#FF69B4", # 热粉
#   "#6A6599", # 紫色
#   "#FFDC91", # 淡黄
#   "#DC0000", # 深红
#   "#708090", # 暗灰蓝
#   "deepskyblue2" , # 天蓝
#   "#AED6F1", # 青绿色
#   "#DEA0FD" # 巧克力棕
# )
#
#
# allcolour <- c(
#   "#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
#   "#6F99AD", "#FFDC91", "#EE4C97", "#3C5488", "#F39B7F",
#   "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
#   "#00A087", "#4DBBD5", "#E64B35", "#374E55", "#DF8F44",
#   "#00A1D5","#E0D4CA","#AB3282","#E95C59","#476D87"
# )

# allcolour <- c(
#   "#C53270", "#57C3F3", "#D6E7A3","lightcoral", "#F1BB72",
#   "tomato3", "#E0D4CA", "#53A85F",  "#00A087","#E59CC4",
#   "#4DBBD5", "#BD956A", "#8C549C",  ,"#5F8670",
#   "#E18727","#DC0000", "tomato4","#7876B1", "#476D87",
#   "#0072B5","#FFD92F")
# allcolour <- c(
#   "#BC3C29", "#0072B5", "#E18727",  "#7E6148","#7876B1",
#   "#BD956A", "#FFDC91", "hotpink3" , "#3C5488", "#F39B7F",
#   "#8491B4", "#91D1C2",  "#6F99AD", "#DC0000","#D6E7A3",
#   "#20854E","#00A087", "#4DBBD5", "#EE4C97" ,"#00A1D5",
#   "violetred4","#FFD92F"
# )

head(merge$seurat_clusters)
head(merge$CellType)


allcolour <- c(
  "#0072B5", "#BC3C29", "#E18727", "#20854E","#00A087", 
  "#6F99AD",  "bisque","lightcoral", "#EE4C97","#3C5488",
  "#D6E7A3", "#7E6148",    "#4DBBD5", "#DC0000", 
  "#7876B1", "#AB3282","#00A1D5","#B19C07","mediumorchid3" ,
  "tomato4","chocolate1","#FFD92F"
)



#cols=alpha(allcolour,0.5)
pdf('umap_name.pdf', width = 8.8, height = 6)
DimPlot(
  merge, 
  reduction = "umap", 
  pt.size = 0.1,  
  label = TRUE,
  raster = FALSE, 
  label.size = 4, 
  cols = allcolour  # 颜色向量，长度 = 细胞类型数
) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=17),
    axis.title.y = element_text(size=17),
    axis.text.x = element_text(size = 15, color="black"),
    axis.text.y = element_text(size = 15, color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="right",
    legend.key.height = unit(0.57,"cm")
  ) +
  guides(
    colour = guide_legend(
      title = "Cell type",
      ncol = 1,                   # 🔹 强制 legend 一列
      override.aes = list(size = 3) # 调整图例点大小
    )
  )

dev.off()


#cols=alpha(allcolour,0.5)
pdf('umap_name-2.pdf', width = 8.8, height = 6)
DimPlot(
  merge, 
  reduction = "umap", 
  pt.size = 2,  
  label = TRUE,
  raster = TRUE, 
  label.size = 4, 
  cols = allcolour  # 颜色向量，长度 = 细胞类型数
) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=17),
    axis.title.y = element_text(size=17),
    axis.text.x = element_text(size = 15, color="black"),
    axis.text.y = element_text(size = 15, color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="right",
    legend.key.height = unit(0.57,"cm")
  ) +
  guides(
    colour = guide_legend(
      title = "Cell type",
      ncol = 1,                   # 🔹 强制 legend 一列
      override.aes = list(size = 3) # 调整图例点大小
    )
  )

dev.off()

#cols=alpha(allcolour,0.5)
pdf('umap_name-3.pdf', width = 8.8, height = 6)
DimPlot(
  merge, 
  reduction = "umap", 
  pt.size = 2,  
  label = F,
  raster = TRUE, 
  label.size = 4, 
  cols = allcolour  # 颜色向量，长度 = 细胞类型数
) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    plot.title = element_text(size=13,hjust = 0.5),
    axis.title.x = element_text(size=17),
    axis.title.y = element_text(size=17),
    axis.text.x = element_text(size = 15, color="black"),
    axis.text.y = element_text(size = 15, color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),       
    panel.background = element_blank(),
    plot.background = element_blank(),   
    legend.background = element_blank(), 
    legend.key        = element_blank(),   
    legend.position="right",
    legend.key.height = unit(0.57,"cm")
  ) +
  guides(
    colour = guide_legend(
      title = "Cell type",
      ncol = 1,                   # 🔹 强制 legend 一列
      override.aes = list(size = 3) # 调整图例点大小
    )
  )

dev.off()

table(merge@meta.data$sample)
# C_1a  C_1b  C_2a  C_2b 
# 9602 27194 21626 27152 

pdf("umap_cluster_type.pdf", width = 20, height = 5)
DimPlot(merge, reduction = "umap",pt.size = 1,  label = F, raster=T, split.by = "sample", label.size = 6 , cols=allcolour) + 
  scale_fill_manual(values = allcolour) +
  xlab('UMAP 1') +
  ylab('UMAP 2')+
  theme_bw()+  
  theme(panel.grid = element_blank(),       
        panel.background = element_blank(),
        plot.background = element_blank() ,   
        legend.background = element_blank(), 
        legend.key        = element_blank(),   
        panel.border = element_rect(fill=NA,color="black", size=1.4, linetype="solid"),
        plot.title = element_text(size=17,hjust = 0.5),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        axis.text.x = element_text(size = 13, color="black"),
        axis.text.y = element_text(size = 13, color="black"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 17),
        legend.position="right",
        legend.key.size = unit(0.46, "cm"),       # 控制每个图例方框的大小
        legend.spacing.y = unit(0.1, "cm"))+
  guides(
    colour = guide_legend(
      title = "Cell type", 
      ncol = 1,            # 🔹 强制一列
      override.aes = list(size = 3)
    )
  )
dev.off()







######################### DimPlot ggplot2
## umap = merge@reductions$umap@cell.embeddings %>% 
umap = merge@reductions$umap@cell.embeddings %>%
  as.data.frame()  %>% 
  cbind(cell = merge@active.ident)
table(umap$cell)

colnames(umap) <- c("umap_1", "umap_2", "cell")
write.csv(umap,"umap.csv")
head(umap)

table <- table(umap$cell) %>% data.frame()
write.csv(table,"cell_level.csv")

##############################################  heatmap #############################################  
# visualize them as a heatmap
# Skin Cell Type Markers - 10 genes each 
# stress response Keratinocytes
stress_response_markers <- c("FOS", "JUN", "JUNB", "ATF3", "EGR1", "DUSP1", "MYC", "NFKBIA", "GADD45B", "KLF4", "CCL2", "CCL27", "IL1R2","HSPA1B" )
# Hyperproliferative Keratinocytes
hyperproliferative_markers <- c("KRT6A","KRT6B","KRT6C", "KRT16", "KRT17")
# Basal Keratinocytes
basal_markers <- c("KRT15","COL17A1","LAMB3","LAMA3","DST")
# Endothelial Cells
blood_endothelial_markers <- c("FLT1", "PLVAP", "CD34","ERG","EMCN", "ESAM")
# Differentiated KCs
# differentiated_KCs_marker <- c("KRT1","KRT2","KRT5","DSC1","CASP14","KLK5")
# Spinous Keratinocytes  (棘层角质形成细胞) 位置: 基底层之上，颗粒层之下
# spinous_markers <- c("PLET1","SPINK5","CLDN1","SBSN","GRHL1","SCEL","LAD1")
# spinous_markers
spinous_KCs <- c("KRT1","KRT2B","KRT77","DSC1","CALML5","SBSN")
# Proliferative Keratinocytes
proliferative_markers  <- c("MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1", "CDK1", "TYMS", "RRM2", "CENPA", "NUSAP1")
# Transitional Keratinocytes 
transitional_markers <- c("IL1R2", "DST")
# Apocrine Sweat Gland Epithelial Cells 大汗腺上皮细胞
apocrine_markers <- c("PIP", "AZGP1", "CRISP3", "CA8", "AQP5", "SCGB2A2" )
# Hair shaft keratinocytes (cortex/cuticle lineage)
hair_shaft_markers <- c("KRT31", "KRT32", "KRT35", "KRT85", "KRT89", "HOXC13",  "DSG4") 
# Sebocytes / Sebaceous Gland Cells 皮脂腺细胞 /Sebaceous Gland Sebocytes (皮脂腺皮脂细胞)
sebocyte_markers <- c("SCD", "THRSP", "ELOVL3", "AWAT1", "AWAT2", "SOAT1", "CIDEA", "PGRMC1", "FASN", "PPARG","KRT7")
# Basal Progenitor_like Keratinocytes 
basal_progenitor_like_markers  <-  c("KRT8","WNT6","LBH","PPARG","IL1R2","DSG2")
# Outer Root Sheath
ors_markers <- c("KRT71", "KRT72", "KRT73", "KRT25", "GATA3", "TCHH", "IVL", "DLX3", "GSDMA", "KRT27")
# Smooth Muscle Cells
smooth_muscle_markers <- c("ACTA2", "MYH11", "TAGLN", "TPM1", "TPM2", "CNN1", "MYLK", "MYL9", "ACTG2", "SMTN")
# Fibroblasts
fibroblast_markers <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "DCN", "LUM", "PDGFRA", "THY1", "FAP")
# Melanocytes
melanocyte_markers <- c("PMEL", "MLANA", "TYRP1", "TYR", "DCT", "MITF", "SOX10", "OCA2", "RAB27A", "GPNMB")
# activated_spinous_KCs
activated_spinous_KCs <- c("ATP6V1C2","UNC93A","KRT6B","KRT17","PLET1","TNFAIP8L3","GRHL1")
# Mast cells
mast_markers <- c("FCER1A", "MS4A2", "HPGDS", "LTC4S")
# T Cells
t_cell_markers <- c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD7", "TRAC", "TRBC1", "LCK", "IL7R")
# Terminally differentiated / cornified keratinocytes
# cornified_markers <- c("KRTDAP", "S100A7", "S100A12", "KLK5", "KLK8", "KLK13", "CASP14", "POU2F3", "DMKN") 
# inflammatory_spinous_KCs
inflammatory_spinous_KCs <- c("S100A7","S100A12","CRABP2","BPIFC","PGLYRP3","FOXQ1")
# Eccrine Sweat Gland Cells 小汗腺细胞
eccrine_markers <- c("AQP5", "ATP1A1", "ATP1B1", "SLC12A2", "EHF")
# Macrophages
macrophage_markers <- c("CD68", "CD163", "AIF1", "FCGR3A", "MARCO", "MSR1", "CD14", "LYZ", "TYROBP", "CTSS")
# Inner Root Sheath
irs_markers <- c("KRT23", "KRT75",  "FOXQ1", "GRHL1", "TCHH", "DSG4", "PADI3", "SERPINB13", "IVL")


markers_genes <- c(blood_endothelial_markers,hyperproliferative_markers,stress_response_markers,basal_markers,
                   hair_shaft_markers,spinous_KCs,transitional_markers, proliferative_markers,apocrine_markers,sebocyte_markers,
                   basal_progenitor_like_markers,smooth_muscle_markers,ors_markers,fibroblast_markers,
                   activated_spinous_KCs,melanocyte_markers,mast_markers,t_cell_markers,inflammatory_spinous_KCs,
                   eccrine_markers,macrophage_markers,irs_markers)

markers_genes <- unique(markers_genes)


allcolour <- c(
  "#0072B5", "#BC3C29", "#E18727", "#20854E","#00A087", 
  "#6F99AD",  "bisque","lightcoral", "#EE4C97","#3C5488",
  "#D6E7A3", "#7E6148",    "#4DBBD5", "#DC0000", 
  "#7876B1", "#AB3282","#00A1D5","#B19C07","mediumorchid3" ,
  "tomato4","chocolate1","#FFD92F"
)


# Build a named vector for group.colors
group_cols <- c(
  "Blood endothelial" ="#0072B5",
  "Hyperproliferative KCs" =  "#BC3C29",
  "Stress response KCs" =  "#E18727", 
  "Basal KCs" =  "#20854E",
  "Hair shaft KCs" = "#00A087", 
  "Spinous KCs" =  "#6F99AD",
  "Transitional KCs" = "bisque", 
  "Proliferative KCs" = "lightcoral",
  "Apocrine sweat gland cells" = "#EE4C97",
  "Sebocytes" = "#3C5488",
  "Basal progenitor-like KCs" = "#D6E7A3",
  "Smooth muscle cells" = "#7E6148",  
  "Outer root sheath KCs" = "#4DBBD5",
  "Fibroblasts" = "#DC0000", 
  "Activated spinous KCs" =  "#7876B1", 
  "Melanocytes" = "#AB3282",
  "Mast cells" = "#00A1D5",
  "T cells" = "#B19C07",
  "Inflammatory spinous KCs" = "mediumorchid3" ,
  "Eccrine sweat gland cells" = "tomato4",
  "Macrophages" =  "chocolate1",
  "Inner root sheath KCs" = "#FFD92F"
)

# create a scale.data slot for the selected genes
merge_1 <- ScaleData(merge, features = unique(c(markers_genes)))
pdf('marker_gene_heatmap.pdf', width = 36, height = 48)
DoHeatmap(merge_1, features = markers_genes, assay = "RNA",group.colors = group_cols)+
  theme(axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
dev.off()




##############################################  VlnPlot #############################################  
# pdf('VlnPlot.pdf', width = 27, height = 54)
# VlnPlot(merge_1, features =markers_genes, stack = TRUE, flip = TRUE, pt.size = 0) +
#   # scale_fill_manual(values = allcolour)+
#   NoLegend() +
#   theme_bw()+  
#   theme(panel.grid=element_blank(), 
#         panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
#         panel.background = element_rect(fill= "white"),
#         plot.title = element_text(size=13,hjust = 0.5),
#         axis.title.x = element_text(size=17),
#         axis.title.y = element_text(size=17),
#         axis.text.x.bottom = element_text(size = 25, angle = 45, hjust = 1),
#         axis.text.y.left = element_blank(), 
#         axis.text.y.right =  element_text(size = 25, angle = 90),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 13),
#         legend.position="none")
# dev.off()


##############################################  DotPlot #############################################  
pdf('DotPlot-2.pdf', width = 48, height = 9)
DotPlot(merge_1, features = markers_genes, assay = "RNA", cols = c("blue", "red")) +
  labs(y="",x="")+
  theme_bw() +
  theme(
    panel.grid        = element_blank(), 
    panel.border      = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid"),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 13, hjust = 0.5),
    axis.title.x      = element_text(size = 17),
    axis.title.y      = element_text(size = 17),
    axis.text.x       = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 13)
  )
dev.off()






##############################################  heatmap #############################################  
# visualize them as a heatmap
# Skin Cell Type Markers - 10 genes each 
stress_response_markers <- c( "MYC", "CCL27", "IL1R2")
# Hyperproliferative Keratinocytes
hyperproliferative_markers <- c("KRT6A", "KRT16", "KRT17")
# Basal Keratinocytes
basal_markers <- c("KRT15","COL17A1","DST")
# Endothelial Cells
blood_endothelial_markers <- c("FLT1", "PLVAP","ERG") #"EMCN", "ESAM"
# Differentiated KCs
# differentiated_KCs_marker <- c("KRT1","KRT2","KRT5","DSC1","CASP14","KLK5")
# Spinous Keratinocytes  (棘层角质形成细胞) 位置: 基底层之上，颗粒层之下
# spinous_markers <- c("PLET1","SPINK5","CLDN1","SBSN","GRHL1","SCEL","LAD1")
# spinous_markers
spinous_KCs <- c("KRT1","KRT2B","KRT77","DSC1")
# Proliferative Keratinocytes
proliferative_markers  <- c("MKI67","TOP2A", "TYMS")
# Transitional Keratinocytes 
transitional_markers <- c("IL1R2", "DST")
# Apocrine Sweat Gland Epithelial Cells 大汗腺上皮细胞
apocrine_markers <- c("PIP", "AZGP1", "SCGB2A2" )
# Hair shaft keratinocytes (cortex/cuticle lineage)
hair_shaft_markers <- c( "KRT35", "KRT89", "HOXC13") 
# Sebocytes / Sebaceous Gland Cells 皮脂腺细胞 /Sebaceous Gland Sebocytes (皮脂腺皮脂细胞)
sebocyte_markers <- c("CIDEA", "ELOVL3", "AWAT1")
# Basal Progenitor_like Keratinocytes 
basal_progenitor_like_markers  <-  c("WNT6","PPARG","DSG2")
# Outer Root Sheath
ors_markers <- c("KRT71", "KRT72", "KRT25")
# Smooth Muscle Cells
smooth_muscle_markers <- c("ACTA2", "MYH11","TPM2")
# Fibroblasts
fibroblast_markers <- c("COL1A1", "FN1",  "DCN")
# Melanocytes
melanocyte_markers <- c("PMEL", "MLANA",  "DCT")
# activated_spinous_KCs
activated_spinous_KCs <- c("ATP6V1C2","UNC93A","PLET1")
# Mast cells
mast_markers <- c("FCER1A", "MS4A2", "HPGDS")
# T Cells
t_cell_markers <- c("CD3D", "CD5", "LCK")
# inflammatory_spinous_KCs
inflammatory_spinous_KCs <- c("S100A7","S100A12","CRABP2")
# Eccrine Sweat Gland Cells 小汗腺细胞
eccrine_markers <- c("ATP1B1", "CLDN8", "EHF")
# Macrophages
macrophage_markers <- c( "AIF1", "TYROBP", "CTSS")
# Inner Root Sheath
irs_markers <- c("KRT75",  "FOXQ1", "SERPINB13")


markers_genes <- c(blood_endothelial_markers,hyperproliferative_markers,stress_response_markers,basal_markers,
                   hair_shaft_markers,spinous_KCs,transitional_markers, proliferative_markers,apocrine_markers,sebocyte_markers,
                   basal_progenitor_like_markers,smooth_muscle_markers,ors_markers,fibroblast_markers,
                   activated_spinous_KCs,melanocyte_markers,mast_markers,t_cell_markers,inflammatory_spinous_KCs,
                   eccrine_markers,macrophage_markers,irs_markers)

markers_genes <- unique(markers_genes)


allcolour <- c(
  "#0072B5", "#BC3C29", "#E18727", "#20854E","#00A087", 
  "#6F99AD",  "bisque","lightcoral", "#EE4C97","#3C5488",
  "#D6E7A3", "#7E6148",    "#4DBBD5", "#DC0000", 
  "#7876B1", "#AB3282","#00A1D5","#B19C07","mediumorchid3" ,
  "tomato4","chocolate1","#FFD92F"
)


# Build a named vector for group.colors
group_cols <- c(
  "Blood endothelial" ="#0072B5",
  "Hyperproliferative KCs" =  "#BC3C29",
  "Stress response KCs" =  "#E18727", 
  "Basal KCs" =  "#20854E",
  "Hair shaft KCs" = "#00A087", 
  "Spinous KCs" =  "#6F99AD",
  "Transitional KCs" = "bisque", 
  "Proliferative KCs" = "lightcoral",
  "Apocrine sweat gland cells" = "#EE4C97",
  "Sebocytes" = "#3C5488",
  "Basal progenitor-like KCs" = "#D6E7A3",
  "Smooth muscle cells" = "#7E6148",  
  "Outer root sheath KCs" = "#4DBBD5",
  "Fibroblasts" = "#DC0000", 
  "Activated spinous KCs" =  "#7876B1", 
  "Melanocytes" = "#AB3282",
  "Mast cells" = "#00A1D5",
  "T cells" = "#B19C07",
  "Inflammatory spinous KCs" = "mediumorchid3" ,
  "Eccrine sweat gland cells" = "tomato4",
  "Macrophages" =  "chocolate1",
  "Inner root sheath KCs" = "#FFD92F"
)


# create a scale.data slot for the selected genes
merge_1 <- ScaleData(merge, features = unique(c(markers_genes)))


pdf('marker_gene_heatmap_short.pdf', width = 36, height = 34)
DoHeatmap(merge_1, features = markers_genes, assay = "RNA",group.colors = group_cols)+
  theme(axis.text.y = element_text(size = 17),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
dev.off()


##############################################  DotPlot #############################################  
pdf('DotPlot_short_right.pdf', width = 25, height = 9)
DotPlot(merge_1, features = markers_genes, assay = "RNA", cols = c("blue", "red")) +
  labs(y="",x="")+
  theme_bw() +
  theme(
    panel.grid        = element_blank(), 
    panel.border      = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid"),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 13, hjust = 0.5),
    axis.title.x      = element_text(size = 17),
    axis.title.y      = element_text(size = 17),
    axis.text.x       = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 13)
  )
dev.off()




pdf('DotPlot_short_top.pdf', width = 23, height = 9)
DotPlot(merge, features = markers_genes, assay = "RNA", cols = c("blue", "red")) +
  labs(y="",x="")+
  theme_bw() +
  theme(
    panel.grid        = element_blank(), 
    panel.border      = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid"),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 13, hjust = 0.5),
    axis.title.x      = element_text(size = 17),
    axis.title.y      = element_text(size = 17),
    axis.text.x       = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 13),
    legend.position = "top"
  )
dev.off()






pdf('DotPlot_short_coord_flip_top.pdf', width = 11, height = 20)
DotPlot(merge, features = markers_genes, assay = "RNA", cols = c("blue", "red")) +
  labs(y="",x="")+
  theme_bw() +
  theme(
    panel.grid        = element_blank(), 
    panel.border      = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid"),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 13, hjust = 0.5),
    axis.title.x      = element_text(size = 17),
    axis.title.y      = element_text(size = 17),
    axis.text.x       = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 13),
    legend.position =   "top"
  )+
  coord_flip()
dev.off()




pdf('DotPlot_short_coord_flip_right.pdf', width = 13, height = 18)
DotPlot(merge, features = markers_genes, assay = "RNA", cols = c("blue", "red")) +
  labs(y="",x="")+
  theme_bw() +
  theme(
    panel.grid        = element_blank(), 
    panel.border      = element_rect(fill = NA, color = "black", size = 1.2, linetype = "solid"),
    panel.background  = element_rect(fill = "white"),
    plot.title        = element_text(size = 13, hjust = 0.5),
    axis.title.x      = element_text(size = 17),
    axis.title.y      = element_text(size = 17),
    axis.text.x       = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 15),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 13),
    legend.position =   "right"
  )+
  coord_flip()
dev.off()





metadata <- merge@meta.data
head(metadata )
write.csv(metadata ,"metadata.csv")


