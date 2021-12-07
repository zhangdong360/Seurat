#使用read_delim函数需要载入readr包
#如未安装，使用
#BiocManager::install('readr')
#安装
library(readr)

# disease -----------------------------------------------------------------


disease_barcodes <- read_delim("disease/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)

View(disease_barcodes)

disease_features <- read_delim("disease/filtered_feature_bc_matrix/features.tsv.gz", 
                       delim = "\t", escape_double = FALSE, 
                       col_names = FALSE, trim_ws = TRUE)

View(disease_features)

#前两行有%注释，可以查看
disease_matrix <- read_delim("disease/filtered_feature_bc_matrix/matrix.mtx.gz", 
                       delim = "\t", escape_double = FALSE, 
                       col_names = FALSE, trim_ws = TRUE)
View(disease_matrix)

#去除前两行注释

disease_matrix <- read.table("D:/sc-data/disease/filtered_feature_bc_matrix/matrix.mtx.gz", quote="\"", comment.char="%")
View(disease_matrix)


# normal ------------------------------------------------------------------
normal_barcodes <- read_delim("normal/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                               delim = "\t", escape_double = FALSE, 
                               col_names = FALSE, trim_ws = TRUE)

View(normal_barcodes)

normal_features <- read_delim("normal/filtered_feature_bc_matrix/features.tsv.gz", 
                               delim = "\t", escape_double = FALSE, 
                               col_names = FALSE, trim_ws = TRUE)

View(normal_features)

#前两行有%注释，可以查看
normal_matrix <- read_delim("normal/filtered_feature_bc_matrix/matrix.mtx.gz", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE)
View(normal_matrix)

#去除前两行注释

normal_matrix <- read.table("normal/filtered_feature_bc_matrix/matrix.mtx.gz", quote="\"", comment.char="%")
View(normal_matrix)


# 后续分析 --------------------------------------------------------------------

BiocManager::install('Seurat')
library(Seurat)
library(dplyr)
library(Matrix)


#读入数据
#disease输入
disease <- Read10X(data.dir = "disease/filtered_feature_bc_matrix/")

#normal输入
disease <- Read10X(data.dir = "normal/filtered_feature_bc_matrix/")
dense.size <- object.size(x = as.matrix(x = disease))
dense.size
#disease 3430246680 bytes
#normal 2846797848 bytes
sparse.size <- object.size(x = disease)
sparse.size
#disease 188236672 bytes
#normal 152268632 bytes

dense.size / sparse.size
#disease 18.2 bytes
#normal 18.7 bytes

# 从原始数据创造Seurat对象
neurons <- CreateSeuratObject(counts = disease, min.cells = 3, min.genes = 200,  project = "10X_neurons")
neurons
#disease
#An object of class Seurat 
#18226 features across 11703 samples within 1 assay 
#Active assay: RNA (18226 features, 0 variable features)

#normal
#An object of class Seurat 
#17620 features across 9711 samples within 1 assay 
#Active assay: RNA (17620 features, 0 variable features)
head(neurons@meta.data)



mito.genes <- grep(pattern = "^MT-", x = rownames(x=neurons), value = TRUE)
percent.mito <- Matrix::colSums(neurons[mito.genes, ]) / Matrix::colSums(neurons)
neurons <- AddMetaData(object = neurons, metadata = percent.mito, col.name = "percent.mito")

#serat2.0 使用以下语句
VlnPlot(object = neurons, features = c( "nGene", "nUMI","percent.mito"), ncol = 3)

#serat3.0 使用以下语句,用它用它
VlnPlot(object = neurons, features = c( "nCount_RNA", "nFeature_RNA","percent.mito"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, 
#but can be used# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

#简书文章里，不同命名的变量
# percent.mt = percent.mito
# pbmc = neurons
plot1 <- FeatureScatter(neurons, feature1 ="nCount_RNA", feature2 ="percent.mito")

plot2 <- FeatureScatter(neurons, feature1 ="nCount_RNA", feature2 ="nFeature_RNA")

CombinePlots(plots = list(plot1, plot2))

#我们过滤具有超过2,500或小于200的独特特征计数的细胞，nFeature_RNA

#我们过滤线粒体计数> 5％的细胞，percent.mt

pbmc <- subset(neurons, subset = nFeature_RNA >200& nFeature_RNA <2500& percent.mito <5)


# 规范化数据

pbmc <- NormalizeData(pbmc, normalization.method ="LogNormalize", scale.factor =10000)


#识别高度可变的 features (FindVariableFeatures)

pbmc <- FindVariableFeatures(pbmc, selection.method ="vst", nfeatures =2000)

# Identify the 10 most highly variable genestop10 <- head(VariableFeatures(pbmc),10)
top10 <- head(VariableFeatures(pbmc),10)
# plot variable features with and without labels

plot1 <- VariableFeaturePlot(pbmc)

plot2 <- LabelPoints(plot = plot1, points = top10, repel =TRUE)

CombinePlots(plots = list(plot1, plot2))

# 缩放数据

all.genes <- rownames(pbmc)
#传输到中间文件pbmc_suofang
pbmc_suofang <- ScaleData(pbmc, features = all.genes)

pbmc <- pbmc_suofang

#执行PCA

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# result
# PC_ 1 
# Positive:  IFI30, CST3, FCN1, LYZ, S100A9, LST1, S100A8, MNDA, CSTA, SERPINA1 
# CD14, MS4A6A, CTSS, SPI1, VCAN, TYMP, FOS, CYBB, FGL2, MPEG1 
# CFD, NCF2, S100A12, CPVL, CD68, FTL, TNFAIP2, LGALS2, HLA-DRA, FTH1 
# Negative:  IFITM1, IL32, CD3D, LDHB, IL7R, PCED1B-AS1, CD52, TCF7, LEF1, EEF1A1 
# CCR7, LTB, CD3G, CD27, NOSIP, CD247, CAMK4, PIK3IP1, LDLRAP1, ITM2A 
# MT-ND4L, RPS18, TRAT1, RPS6, FLT3LG, MYC, RPS12, MAL, ABLIM1, PRKCQ-AS1 
# PC_ 2 
# Positive:  NKG7, CST7, PRF1, FGFBP2, GNLY, GZMB, GZMH, GZMA, CTSW, KLRD1 
# CCL5, FCGR3A, SPON2, HOPX, CX3CR1, CLIC3, KLRC2, ADGRG1, KLRF1, PLEK 
# SRGN, S1PR5, PRSS23, CCL4, MYOM2, TTC38, TYROBP, FCRL6, CD247, APMAP 
# Negative:  LTB, RPS8, RPL32, CCR7, RPL13, TPT1, RPS12, RPLP1, LEF1, SELL 
# TCF7, RPLP0, IL7R, RPS18, EEF1A1, RPS6, LDHB, CD27, MYC, CAMK4 
# NOSIP, ACTN1, LDLRAP1, MAL, RCAN3, PIK3IP1, TRAT1, COTL1, NCF1, TRABD2A 
# PC_ 3 
# Positive:  MS4A1, CD79A, BANK1, CD79B, TCL1A, HLA-DQA1, NIBAN3, IGHM, FCRLA, SPIB 
# LINC00926, HLA-DOB, HLA-DQB1, VPREB3, CD22, FCER2, CD19, RALGPS2, BLNK, CD74 
# FCRL1, HLA-DPB1, CD24, TSPAN13, HLA-DRA, SWAP70, HVCN1, CD72, HLA-DRB1, IGHD 
# Negative:  IL7R, TCF7, LEF1, CD3D, LDHB, AIF1, TPT1, NOSIP, CD3G, ACTN1 
# CAMK4, S100A6, IL32, SAMHD1, VIM, TRAT1, IFITM1, LDLRAP1, MAL, S100A8 
# S100A10, S100A9, FCN1, RPLP1, CD27, RPS12, PIK3IP1, LYZ, CD14, PRKCQ-AS1 
# PC_ 4 
# Positive:  LILRA4, PLD4, IL3RA, SERPINF1, CLEC4C, LRRC26, PPP1R14B, SCT, MAP1A, UGCG 
# PTCRA, ITM2C, APP, TPM2, EPHB1, IRF7, TRAF4, PTPRS, CCDC50, LILRB4 
# IRF8, PPM1J, SMPD3, LAMP5, LINC00996, TLR9, P2RY14, TNFRSF21, FCER1A, TCF4 
# Negative:  IGLV2-14, IGHV7-4-1, IGHV3-30, IGKV3-20, S100A8, IGKV3-15, IGHV4-39, HBB, S100A9, IGHV3-15 
# IGKV4-1, LYZ, IGHV4-59, IGKV2-30, IGLV3-19, IGKV1-16, VCAN, MS4A1, CD79A, GNLY 
# IGHV3-74, MALAT1, S100A12, FCER2, VPREB3, LINC00926, MNDA, CD22, FCN1, IGHD 
# PC_ 5 
# Positive:  PPBP, TUBB1, JCHAIN, CAVIN2, GNG11, SPARC, PF4, F13A1, NRGN, IGHV7-4-1 
# MPIG6B, IGLV2-14, CLU, IGKV3-20, IGHV3-30, GP9, IGKV3-15, PTGS1, MYL9, PTCRA 
# GP1BB, RGS18, HBB, IGHV4-39, IGHV3-15, ITGA2B, TSC22D1, PGRMC1, CMTM5, IGKV2-30 
# Negative:  EEF1A1, CD52, RPS18, RPL13, RPS6, RPS12, RPL32, JUNB, RPLP1, RPS8 
# CNN2, CD37, RPLP0, IL32, CD3D, LTB, CRIP1, DUSP1, IFITM2, TPT1 
# DUSP2, IFITM1, GZMK, CD3G, S100A10, EZR, ANXA1, ACTG1, KLRG1, CYBA 

#提供定义可视化细胞和功能的几种有用的方式PCA，包括VizDimReduction，DimPlot，和DimHeatmap
print(pbmc[["pca"]], dims =1:5, nfeatures =5)
# 1 
# Positive:  IFI30, CST3, FCN1, LYZ, S100A9 
# Negative:  IFITM1, IL32, CD3D, LDHB, IL7R 
# PC_ 2 
# Positive:  NKG7, CST7, PRF1, FGFBP2, GNLY 
# Negative:  LTB, RPS8, RPL32, CCR7, RPL13 
# PC_ 3 
# Positive:  MS4A1, CD79A, BANK1, CD79B, TCL1A 
# Negative:  IL7R, TCF7, LEF1, CD3D, LDHB 
# PC_ 4 
# Positive:  LILRA4, PLD4, IL3RA, SERPINF1, CLEC4C 
# Negative:  IGLV2-14, IGHV7-4-1, IGHV3-30, IGKV3-20, S100A8 
# PC_ 5 
# Positive:  PPBP, TUBB1, JCHAIN, CAVIN2, GNG11 
# Negative:  EEF1A1, CD52, RPS18, RPL13, RPS6

## VizDimLoadings方法
VizDimLoadings(pbmc, dims =1:2, reduction ="pca")

## DimPlot方法

DimPlot(pbmc, reduction ="pca")

## DimHeatmap方法
DimHeatmap(pbmc, dims =1, cells =500, balanced =TRUE)
### 输出PC1至PC15（个人推测）
DimHeatmap(pbmc, dims =1:15, cells =500, balanced =TRUE)

# 确定数据集的“维度”

# NOTE:This process can take a long time for big datasets, comment out for expediency. 
# More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
# 翻译：#注意:对于大数据集，这个过程可能需要很长时间，为了方便起见，请注释掉。 
# 可以使用更近似的技术(如ElbowPlot()中实现的技术)来减少计算时间  

pbmc <- JackStraw(pbmc, num.replicate =100)

pbmc <- ScoreJackStraw(pbmc, dims =1:20)

JackStrawPlot(pbmc, dims =1:20)

# 另一种确定维度的方法，即方差（ElbowPlot函数）
ElbowPlot(pbmc)

# 细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

# 非线性降维分析
pbmc_PCA <- RunUMAP(pbmc, dims =1:20)

DimPlot(pbmc_PCA, reduction ="umap")


## 添加细胞类群xiba的标签
DimPlot(pbmc_PCA, reduction = "umap",label = TRUE)
#或者使用
#LabelClusters(DimPlot(pbmc_PCA, reduction = "umap"),id = 'ident')

saveRDS(pbmc_PCA, file ="result/normal/pbmc_tutorial.rds")

## 寻找差异表达的特征（聚类生物标志物）

cluster1.markers <- FindMarkers(pbmc_PCA, ident.1 =1, min.pct =0.25)

head(cluster1.markers, n =5)
# 输出结果至csv文件
write.csv(cluster1.markers,file = "result/normal/cluster1.markers.csv")

# 找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞  
pbmc.markers <- FindAllMarkers(pbmc_PCA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(pbmc.markers,file = 'pbmc.markers.csv')



VlnPlot(pbmc, features = c("S100A9", "HBB"))
FeaturePlot(pbmc_PCA, features = c("S100A9","HBB"))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc_PCA, features = top10$gene) + NoLegend()

# 确定细胞聚类簇的种类
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

