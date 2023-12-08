library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "~/Desktop/Eye/GSE148063")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(pbmc)
plot1

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)

#pbmc <- RunUMAP(pbmc, dims = 1:15)
#DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "tsne")

# write(pbmc@meta.data, "stage.csv")
# Add one empty column and leave single row indicating stage
stage <- read.csv("stage.csv", header = FALSE)
CellsMeta = pbmc@meta.data
head(CellsMeta)
CellsMeta["stage"] <- stage
head(CellsMeta)
pbmc <- AddMetaData(pbmc, CellsMeta)
pbmc
pbmc@meta.data

saveRDS(pbmc, file = "../pbmc_QC.rds")
# pbmc <- readRDS("pbmc_QC.rds")

# png("tsne.png")
# DimPlot(pbmc, reduction = "tsne")
# dev.off()    
# 
# png("tsne_stage.png")
# DimPlot(pbmc, reduction = "tsne", group.by = "stage")
# dev.off()
# 
# png("Notch_tsne.png")
# FeaturePlot(pbmc, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
# dev.off()
# 
# png("Delta_tsne.png")
# FeaturePlot(pbmc, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
# dev.off()
# 
# png("target_tsne.png")
# FeaturePlot(pbmc, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
# dev.off()
# 
# png("mitosis_tsne.png")
# FeaturePlot(pbmc, features = c("Mki67"))
# dev.off()
# 
# # optic vesicle markers "Pax6", "Six3", "Rax", "Lhx2"
# png("RPC_tsne.png")
# FeaturePlot(pbmc, features = c("Pax6", "Six3", "Rax", "Lhx2"))
# dev.off()
# 
# # mesenchymal markers "Twist1", "Col3a1", "Alx3", "Prrx2"
# png("mes_tsne.png")
# FeaturePlot(pbmc, features = c("Twist1", "Col3a1", "Alx3", "Prrx2"))
# dev.off()
# 
# # surface ectoderm markers "Krt18", "Krt8", "Dlx5", "Dlx6"
# png("se_tsne.png")
# FeaturePlot(pbmc, features = c("Krt18"))
# dev.off()
# 
# # endothelium markers "Pecam1", "Cd34", "Cldn5", "Cdh5"
# png("endo_tsne.png")
# FeaturePlot(pbmc, features = c("Pecam1"))
# dev.off()
# 
# # erythrocyte markers "Hbb-bs", "Hbb-bt", "Alas2", "Gypa"
# png("eryth_tsne.png")
# FeaturePlot(pbmc, features = c("Hbb-bs"))
# dev.off()
# 
# # RPC Lhx2, mes Twist1, se Krt18, endo Pecam1, eryth Hbb-bs
# cd_genes <- c("Lhx2", "Twist1", "Krt18", "Pecam1", "Hbb-bs")
# 
# png("dot_marker.png")
# DotPlot(object = pbmc, features = cd_genes)
# dev.off()
# 
# head(pbmc)
# 
new.cluster.ids <- c("RPC1", "RPC2", "RPC3","mes1", "mes2", "mes3", "RPC4", "mes4", "ectoderm", "RPC5", "endothelium", "RPC6", "RPC7", "erythrocyte", "no ident 1", "no ident 2", "no ident 3") 
names(new.cluster.ids) <- levels(pbmc)
pbmc_whole <- RenameIdents(pbmc, new.cluster.ids)
# 
# png("annotation.png")
# DimPlot(pbmc_whole, reduction = "tsne", label = TRUE)
# dev.off()
# 
# png("tsne.png")
# DimPlot(pbmc_whole, reduction = "tsne", label = FALSE)
# dev.off() 

DimPlot(pbmc_whole, reduction = "tsne", label = TRUE, split.by = "stage")

DotPlot(object = pbmc_whole, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = pbmc_whole, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = pbmc_whole, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

FeaturePlot(pbmc_whole, features = "Jag1", cols = c("gray", "red"))
FeaturePlot(pbmc_whole, features = "Dll1", cols = c("gray", "red"))
FeaturePlot(pbmc_whole, features = "Notch1", cols = c("gray", "red"))
FeaturePlot(pbmc_whole, features = "Notch3", cols = c("gray", "red"))
FeaturePlot(pbmc_whole, features = "Hes1", cols = c("gray", "red"))

FeaturePlot(pbmc_whole, features = "Jag1", cols = c("gray", "red"), split.by = "stage", label = TRUE)
FeaturePlot(pbmc_whole, features = "Dll1", cols = c("gray", "red"), split.by = "stage", label = TRUE)
FeaturePlot(pbmc_whole, features = "Notch1", cols = c("gray", "red"), split.by = "stage", label = TRUE)
FeaturePlot(pbmc_whole, features = "Notch3", cols = c("gray", "red"), split.by = "stage", label = TRUE)
FeaturePlot(pbmc_whole, features = "Hes1", cols = c("gray", "red"), split.by = "stage", label = TRUE)

VlnPlot(pbmc_whole, features = "Jag1")
VlnPlot(pbmc_whole, features = "Dll1")
VlnPlot(pbmc_whole, features = "Notch1")
VlnPlot(pbmc_whole, features = "Notch3")

pbmc_whole_RPC <- subset(pbmc_whole, idents = c("RPC1", "RPC2", "RPC3","RPC4", "RPC5", "RPC6", "RPC7"))
# pbmc_whole_RPC 4508 samples

# RPC1.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC1")
# RPC2.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC2")
# RPC3.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC3")
# RPC4.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC4")
# RPC5.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC5")
# RPC6.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC6")
# RPC7.markers <- FindMarkers(pbmc_whole_RPC, ident.1 = "RPC7")
# write.csv(RPC1.markers, "RPC1_markers.csv")
# write.csv(RPC2.markers, "RPC2_markers.csv")
# write.csv(RPC3.markers, "RPC3_markers.csv")
# write.csv(RPC4.markers, "RPC4_markers.csv")
# write.csv(RPC5.markers, "RPC5_markers.csv")
# write.csv(RPC6.markers, "RPC6_markers.csv")
# write.csv(RPC7.markers, "RPC7_markers.csv")

# Dmrta2 is involved in NPC maintenance PMID: 28655839
# Wnt8b supresses neuroretina specification PMID: 20890044
VlnPlot(pbmc_whole_RPC, features = "Wnt8b") 
DotPlot(pbmc_whole_RPC, features = "Wnt8b")
DotPlot(pbmc_whole_RPC, features = c("Dmrta2", "Wnt8b", "Sox2") )

VlnPlot(pbmc_whole_RPC, features = "Jag1")
VlnPlot(pbmc_whole_RPC, features = "Dll1")
VlnPlot(pbmc_whole_RPC, features = "Notch1")
VlnPlot(pbmc_whole_RPC, features = "Notch3")

DotPlot(pbmc_whole, features = c("Pax6", "Six3", "Rax", "Lhx2"))
DotPlot(pbmc_whole_RPC, features = c("Olig2", "Ascl1", "Neurog2", "Cdh6", "Otx2", "Vsx1", "Oc1"))

pbmc_whole_RPC_Jag1 <- subset(x = pbmc_whole_RPC, subset =  Jag1 > 0)
# pbmc_whole_RPC_Jag1 905 samples (20.07 %)

pbmc_whole_RPC_Dll1 <- subset(x = pbmc_whole_RPC, subset =  Dll1 > 0)
# pbmc_whole_RPC_Dll1 563 samples (12.48 %)

pbmc_whole_RPC_Notch3 <- subset(x = pbmc_whole_RPC, subset =  Notch3 > 0)
# pbmc_whole_RPC_Notch3 1978 samples (43.87 %)

pbmc_whole_RPC_Notch1 <- subset(x = pbmc_whole_RPC, subset =  Notch1 > 0)
# pbmc_whole_RPC_Notch1 1567 samples (34.76 %)

pbmc_whole_RPC_Hes1 <- subset(x = pbmc_whole_RPC, subset =  Hes1 > 0)
# pbmc_whole_RPC_Hes1 4157 samples (92.21 %)

FeatureScatter(pbmc_whole_RPC, "Notch3", "Hes1")
pbmc_whole_RPC_Notch3_Hes1 <- subset(x = pbmc_whole_RPC, subset = Notch3 > 0 & Hes1 > 0)
# pbmc_whole_RPC_Notch3_Hes1 1896 samples (42.05 %)

FeatureScatter(pbmc_whole_RPC, "Notch1", "Hes1")
pbmc_whole_RPC_Notch1_Hes1 <- subset(x = pbmc_whole_RPC, subset = Notch1 > 0 & Hes1 > 0)
# pbmc_whole_RPC_Notch1_Hes1 1494 samples (33.14 %)

FeatureScatter(pbmc_whole_RPC, "Notch1", "Jag1")
pbmc_whole_RPC_Notch1_Jag1 <- subset(x = pbmc_whole_RPC, subset = Notch1 > 0 & Jag1 > 0)
# pbmc_whole_RPC_Notch1_Jag1 370 samples (8.20 %) (rest 91.79 %)

FeatureScatter(pbmc_whole_RPC, "Notch1", "Dll1")
pbmc_whole_RPC_Notch1_Dll1 <- subset(x = pbmc_whole_RPC, subset = Notch1 > 0 & Dll1 > 0)
# pbmc_whole_RPC_Notch1_Dll1 249 samples (5.52 %) (rest 94.47 %)

FeatureScatter(pbmc_whole_RPC, "Notch3", "Jag1")
pbmc_whole_RPC_Notch3_Jag1 <- subset(x = pbmc_whole_RPC, subset = Notch3 > 0 & Jag1 > 0)
# pbmc_whole_RPC_Notch3_Jag1 494 samples (10.95 %) (rest 89.04 %)

FeatureScatter(pbmc_whole_RPC, "Notch3", "Dll1")
pbmc_whole_RPC_Notch3_Dll1 <- subset(x = pbmc_whole_RPC, subset = Notch3 > 0 & Dll1 > 0)
# pbmc_whole_RPC_Notch3_Dll1 252 samples (5.59 %) (rest 94.40 %)

pbmc_whole_RPC_Jag1_Wnt8b <- subset(x = pbmc_whole_RPC, subset = Jag1 > 0 & Wnt8b > 0)
# pbmc_whole_RPC_Jag1_Wnt8b 183 samples

#################################

somite_12 <- subset(x = pbmc, cells = WhichCells(pbmc, expression = stage =="somite_12"))
somite_16 <- subset(x = pbmc, cells = WhichCells(pbmc, expression = stage =="somite_16"))
somite_20 <- subset(x = pbmc, cells = WhichCells(pbmc, expression = stage =="somite_20"))
somite_24 <- subset(x = pbmc, cells = WhichCells(pbmc, expression = stage =="somite_24"))
somite_26 <- subset(x = pbmc, cells = WhichCells(pbmc, expression = stage =="somite_26"))

DimPlot(somite_12, reduction = "tsne")
new.cluster.ids.12 <- c(
                        "RPC2", # 1
                        "mes2", # 4
                        "RPC4",  # 6
                        "mes4", # 7
                        "ectoderm", # 8
                        "RPC5", # 9 
                        "endothelium", # 10
                        "RPC6", # 11
                        "RPC7", # 12
                        "erythrocyte", # 13
                        "no ident 1", # 14
                        "no ident 3"  # 16
) 
names(new.cluster.ids.12) <- levels(somite_12)
somite_12 <- RenameIdents(somite_12, new.cluster.ids.12)
DimPlot(somite_12, reduction = "tsne", label = TRUE)

DimPlot(somite_16, reduction = "tsne")
new.cluster.ids.16 <- c("RPC1", # 0
                        "RPC2", # 1
                        "RPC3", # 2
                        "mes1", # 3
                        "mes2", # 4
                        "RPC4",  # 6
                        "mes4", # 7
                        "ectoderm", # 8
                        "RPC5", # 9 
                        "endothelium", # 10
                        "RPC6", # 11
                        "RPC7", # 12
                        "erythrocyte", # 13
                        "no ident 1", # 14
                        "no ident 3"  # 16
) 
names(new.cluster.ids.16) <- levels(somite_16)
somite_16 <- RenameIdents(somite_16, new.cluster.ids.16)
DimPlot(somite_16, reduction = "tsne", label = TRUE)

DimPlot(somite_20, reduction = "tsne")
new.cluster.ids.20 <- c("RPC1", # 0
                        "RPC2", # 1
                        "RPC3", # 2
                        "mes1", # 3
                        "mes2", # 4
                        "RPC4",  # 6
                        "mes4", # 7
                        "ectoderm", # 8
                        "RPC5", # 9 
                        "endothelium", # 10
                        "RPC6", # 11
                        "RPC7", # 12
                        "erythrocyte", # 13
                        "no ident 1", # 14
                        "no ident 2", # 15 
                        "no ident 3"  # 16
) 
names(new.cluster.ids.20) <- levels(somite_20)
somite_20 <- RenameIdents(somite_20, new.cluster.ids.20)
DimPlot(somite_20, reduction = "tsne", label = TRUE)

DimPlot(somite_24, reduction = "tsne")
new.cluster.ids.24 <- c("RPC1", # 0
                        "RPC2", # 1
                        "RPC3", # 2
                        "mes1", # 3
                        "mes2", # 4
                        "mes3", # 5
                        "RPC4",  # 6
                        "mes4", # 7
                        "ectoderm", # 8
                        "RPC5", # 9 
                        "endothelium", # 10
                        "RPC6", # 11
                        "RPC7", # 12
                        "erythrocyte", # 13
                        "no ident 1", # 14
                        "no ident 2", # 15 
                        "no ident 3"  # 16
) 
names(new.cluster.ids.24) <- levels(somite_24)
somite_24 <- RenameIdents(somite_24, new.cluster.ids.24)
DimPlot(somite_24, reduction = "tsne", label = TRUE)

DimPlot(somite_26, reduction = "tsne")
new.cluster.ids.26 <- c("RPC1", # 0
                        "RPC2", # 1
                        "RPC3", # 2
                        "mes1", # 3
                        "mes3", # 5
                        "RPC4",  # 6
                        "mes4", # 7
                        "ectoderm", # 8
                        "RPC5", # 9 
                        "endothelium", # 10
                        "RPC6", # 11
                        "RPC7", # 12
                        "erythrocyte", # 13
                        "no ident 1", # 14
                        "no ident 2", # 15 
                        "no ident 3"  # 16
) 
names(new.cluster.ids.26) <- levels(somite_26)
somite_26 <- RenameIdents(somite_26, new.cluster.ids.26)
DimPlot(somite_26, reduction = "tsne", label = TRUE)

# pbmc_somite_12_RPC1 <- subset(somite_12, idents = "RPC1") # n = 0
pbmc_somite_12_RPC2 <- subset(somite_12, idents = "RPC2") # n = 544
# pbmc_somite_12_RPC3 <- subset(somite_12, idents = "RPC3") # n = 0
pbmc_somite_12_RPC4 <- subset(somite_12, idents = "RPC4") # n = 51
pbmc_somite_12_RPC5 <- subset(somite_12, idents = "RPC5") # n = 17
pbmc_somite_12_RPC6 <- subset(somite_12, idents = "RPC6") # n = 29
pbmc_somite_12_RPC7 <- subset(somite_12, idents = "RPC7") # n = 3

pbmc_somite_16_RPC1 <- subset(somite_16, idents = "RPC1") # n = 15
pbmc_somite_16_RPC2 <- subset(somite_16, idents = "RPC2") # n = 441
pbmc_somite_16_RPC3 <- subset(somite_16, idents = "RPC3") # n = 12
pbmc_somite_16_RPC4 <- subset(somite_16, idents = "RPC4") # n = 155
pbmc_somite_16_RPC5 <- subset(somite_16, idents = "RPC5") # n = 153
pbmc_somite_16_RPC6 <- subset(somite_16, idents = "RPC6") # n = 92
pbmc_somite_16_RPC7 <- subset(somite_16, idents = "RPC7") # n = 26

pbmc_somite_20_RPC1 <- subset(somite_20, idents = "RPC1") # n = 667
pbmc_somite_20_RPC2 <- subset(somite_20, idents = "RPC2") # n = 8
pbmc_somite_20_RPC3 <- subset(somite_20, idents = "RPC3") # n = 1
pbmc_somite_20_RPC4 <- subset(somite_20, idents = "RPC4") # n = 146
pbmc_somite_20_RPC5 <- subset(somite_20, idents = "RPC5") # n = 151
pbmc_somite_20_RPC6 <- subset(somite_20, idents = "RPC6") # n = 95
pbmc_somite_20_RPC7 <- subset(somite_20, idents = "RPC7") # n = 110

pbmc_somite_24_RPC1 <- subset(somite_24, idents = "RPC1") # n = 328
pbmc_somite_24_RPC2 <- subset(somite_24, idents = "RPC2") # n = 1
pbmc_somite_24_RPC3 <- subset(somite_24, idents = "RPC3") # n = 101
pbmc_somite_24_RPC4 <- subset(somite_24, idents = "RPC4") # n = 130
pbmc_somite_24_RPC5 <- subset(somite_24, idents = "RPC5") # n = 82
pbmc_somite_24_RPC6 <- subset(somite_24, idents = "RPC6") # n = 91
pbmc_somite_24_RPC7 <- subset(somite_24, idents = "RPC7") # n = 96

pbmc_somite_26_RPC1 <- subset(somite_26, idents = "RPC1") # n = 32
pbmc_somite_26_RPC2 <- subset(somite_26, idents = "RPC2") # n = 5
pbmc_somite_26_RPC3 <- subset(somite_26, idents = "RPC3") # n = 802
pbmc_somite_26_RPC4 <- subset(somite_26, idents = "RPC4") # n = 2
pbmc_somite_26_RPC5 <- subset(somite_26, idents = "RPC5") # n = 5
pbmc_somite_26_RPC6 <- subset(somite_26, idents = "RPC6") # n = 68
pbmc_somite_26_RPC7 <- subset(somite_26, idents = "RPC7") # n = 49

#################################

pbmc_somite_12_RPC2_Wnt8b <- subset(x = pbmc_somite_12_RPC2, 
                                   subset = Wnt8b > 0)
pbmc_somite_12_RPC2_Six3 <- subset(x = pbmc_somite_12_RPC2, 
                                   subset = Six3 > 0)
pbmc_somite_12_RPC2_Jag1 <- subset(x = pbmc_somite_12_RPC2, 
                                        subset = Jag1 > 0)
pbmc_somite_12_RPC2_Dll1 <- subset(x = pbmc_somite_12_RPC2, 
                                  subset = Dll1 > 0)
pbmc_somite_12_RPC2_Notch1 <- subset(x = pbmc_somite_12_RPC2, 
                                  subset = Notch1 > 0)
pbmc_somite_12_RPC2_Notch3 <- subset(x = pbmc_somite_12_RPC2, 
                                    subset = Notch3 > 0)

pbmc_somite_12_RPC2_Jag1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                        subset = Jag1 > 0)
pbmc_somite_12_RPC2_Dll1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                        subset = Dll1 > 0)
pbmc_somite_12_RPC2_Notch1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                        subset = Notch1 > 0)
pbmc_somite_12_RPC2_Notch3_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                          subset = Notch3 > 0)
pbmc_somite_12_RPC2_Six3_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                          subset = Six3 > 0)
pbmc_somite_12_RPC2_Jag1_Dll1 <- subset(x = pbmc_somite_12_RPC2, 
                                        subset = Jag1 > 0 & Dll1 > 0)

pbmc_somite_12_RPC2_Jag1_Notch1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                        subset = Jag1 > 0 & Notch1 > 0)
pbmc_somite_12_RPC2_Jag1_Notch3_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                               subset = Jag1 > 0 & Notch3 > 0)
pbmc_somite_12_RPC2_Dll1_Notch1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                               subset = Dll1 > 0 & Notch1 > 0)
pbmc_somite_12_RPC2_Dll1_Notch3_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                               subset = Dll1 > 0 & Notch3 > 0)
pbmc_somite_12_RPC2_Jag1_Dll1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Wnt8b, 
                                               subset = Jag1 > 0 & Dll1 > 0)

pbmc_somite_12_RPC2_Jag1_or_Dll1 <- subset(x = pbmc_somite_12_RPC2, 
                                          subset = Jag1 > 0 | Dll1 > 0)

pbmc_somite_12_RPC2_Jag1_or_Dll1_Wnt8b <- subset(x = pbmc_somite_12_RPC2_Jag1_or_Dll1, 
                                          subset = Wnt8b > 0)

# pbmc_somite_12_RPC2 n = 544

# pbmc_somite_12_RPC2_Wnt8b n = 235
# pbmc_somite_12_RPC2_Six3 n = 392
# pbmc_somite_12_RPC2_Jag1 n = 163
# pbmc_somite_12_RPC2_Dll1 n = 160
# pbmc_somite_12_RPC2_Notch1 n = 108
# pbmc_somite_12_RPC2_Notch3 n = 200

# pbmc_somite_12_RPC2_Jag1_Wnt8b n = 83
# pbmc_somite_12_RPC2_Dll1_Wnt8b n = 94
# pbmc_somite_12_RPC2_Notch1_Wnt8b n = 51
# pbmc_somite_12_RPC2_Notch3_Wnt8b n = 85
# pbmc_somite_12_RPC2_Six3_Wnt8b n = 142
# pbmc_somite_12_RPC2_Jag1_Dll1 n = 59

# pbmc_somite_12_RPC2_Jag1_Notch1_Wnt8b n = 30
# pbmc_somite_12_RPC2_Jag1_Notch3_Wnt8b n = 36
# pbmc_somite_12_RPC2_Dll1_Notch1_Wnt8b n = 23
# pbmc_somite_12_RPC2_Dll1_Notch3_Wnt8b n = 35
# pbmc_somite_12_RPC2_Jag1_Dll1_Wnt8b n = 39

# pbmc_somite_12_RPC2_Jag1_or_Dll1 n = 264
# pbmc_somite_12_RPC2_Jag1_or_Dll1_Wnt8b n = 138

#################################

pbmc_somite_16_RPC2_Wnt8b <- subset(x = pbmc_somite_16_RPC2, 
                                   subset = Wnt8b > 0)

pbmc_somite_16_RPC2_Jag1 <- subset(x = pbmc_somite_16_RPC2, 
                                  subset = Jag1 > 0)
pbmc_somite_16_RPC2_Dll1 <- subset(x = pbmc_somite_16_RPC2, 
                                  subset = Dll1 > 0)
pbmc_somite_16_RPC2_Jag1_Wnt8b <- subset(x = pbmc_somite_16_RPC2, 
                                  subset = Jag1 > 0 & Wnt8b > 0)
pbmc_somite_16_RPC2_Dll1_Wnt8b <- subset(x = pbmc_somite_16_RPC2, 
                                  subset = Dll1 > 0 & Wnt8b > 0)
pbmc_somite_16_RPC2_Jag1_or_Dll1 <- subset(x = pbmc_somite_16_RPC2, 
                                          subset = Jag1 > 0 | Dll1 > 0)
pbmc_somite_16_RPC2_Jag1_or_Dll1_Wnt8b <- subset(x = pbmc_somite_16_RPC2_Jag1_or_Dll1, 
                                          subset = Wnt8b > 0)

# pbmc_somite_16_RPC2 n = 441
# pbmc_somite_16_RPC2_Wnt8b n = 65
# pbmc_somite_16_RPC2_Jag1 n = 106
# pbmc_somite_16_RPC2_Dll1 n = 84
# pbmc_somite_16_RPC2_Jag1_Wnt8b n = 17
# pbmc_somite_16_RPC2_Dll1_Wnt8b n = 7
# pbmc_somite_16_RPC2_Jag1_or_Dll1 n = 171
# pbmc_somite_16_RPC2_Jag1_or_Dll1_Wnt8b n = 24
#################################

FeaturePlot(somite_12, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = somite_12, features = c("Notch1", "Notch2", "Notch3", "Notch4"))

FeaturePlot(somite_12, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = somite_12, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))

FeaturePlot(somite_12, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
DotPlot(object = somite_12, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

somite_12_RPC <- subset(somite_12, idents = c("RPC2", "RPC4", "RPC5", "RPC6", "RPC7"))
FeatureScatter(somite_12_RPC, "Notch3", "Hes1")
FeatureScatter(somite_12_RPC, "Notch1", "Hes1")
FeatureScatter(somite_12_RPC, "Jag1", "Hes1")
FeatureScatter(somite_12_RPC, "Dll1", "Hes1")
FeatureScatter(somite_12_RPC, "Notch3", "Dll1")

somite_12_RPC_Hes1 <- subset(somite_12_RPC, subset = Hes1 > 0)
FeatureScatter(somite_12_RPC_Hes1, "Notch3", "Notch1")
FeatureScatter(somite_12_RPC_Hes1, "Notch3", "Dll1")

somite_12_RPC_Notch3 <- subset(somite_12_RPC, subset = Notch3 > 0)
FeatureScatter(somite_12_RPC_Notch3, "Dll1", "Hes1")

#################################

FeaturePlot(somite_16, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = somite_16, features = c("Notch1", "Notch2", "Notch3", "Notch4"))

FeaturePlot(somite_16, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = somite_16, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))

FeaturePlot(somite_16, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
DotPlot(object = somite_16, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

somite_16_RPC <- subset(somite_16, idents = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7"))
FeatureScatter(somite_16_RPC, "Notch3", "Hes1")
FeatureScatter(somite_16_RPC, "Notch1", "Hes1")

#################################

FeaturePlot(somite_20, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = somite_20, features = c("Notch1", "Notch2", "Notch3", "Notch4"))

FeaturePlot(somite_20, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = somite_20, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))

FeaturePlot(somite_20, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
DotPlot(object = somite_20, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

somite_20_RPC <- subset(somite_20, idents = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7"))
FeatureScatter(somite_20_RPC, "Notch3", "Hes1")
FeatureScatter(somite_20_RPC, "Notch1", "Hes1")

#################################

FeaturePlot(somite_24, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = somite_24, features = c("Notch1", "Notch2", "Notch3", "Notch4"))

FeaturePlot(somite_24, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = somite_24, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))

FeaturePlot(somite_24, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
DotPlot(object = somite_24, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

somite_24_RPC <- subset(somite_24, idents = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7"))
FeatureScatter(somite_24_RPC, "Notch3", "Hes1")
FeatureScatter(somite_24_RPC, "Notch1", "Hes1")

#################################

FeaturePlot(somite_26, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
DotPlot(object = somite_26, features = c("Notch1", "Notch2", "Notch3", "Notch4"))

FeaturePlot(somite_26, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))
DotPlot(object = somite_26, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"))

FeaturePlot(somite_26, features = c("Hes1", "Hes5", "Hey1", "Hey2"))
DotPlot(object = somite_26, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

somite_26_RPC <- subset(somite_26, idents = c("RPC1", "RPC2", "RPC3", "RPC4", "RPC5", "RPC6", "RPC7"))
FeatureScatter(somite_26_RPC, "Notch3", "Hes1")
FeatureScatter(somite_26_RPC, "Notch1", "Hes1")

#################################

new.cluster.ids_2 <- c("RPC", "RPC", "RPC","mes", "mes", "mes", "RPC", "mes", "endoderm", "RPC", "endothelium", "RPC", "RPC", "erythrocyte", "no ident", "no ident", "no ident") 
names(new.cluster.ids_2) <- levels(pbmc)
pbmc_integrate <- RenameIdents(pbmc, new.cluster.ids_2)
png("annotation_RPC.png")
DimPlot(pbmc_integrate, reduction = "tsne", label = TRUE)
dev.off()

pbmc_RPC <- subset(pbmc_integrate, idents = c("RPC"))
head(pbmc_RPC)

png("RPC_stage.png")
DimPlot(pbmc_RPC, reduction = "tsne", group.by = "stage")
dev.off()

png("dot_RPC_Notch.png")
DotPlot(object = pbmc_RPC, features = c("Notch1", "Notch2", "Notch3", "Notch4"), group.by = "stage")
dev.off()

png("dot_RPC_Delta.png")
DotPlot(object = pbmc_RPC, features = c("Jag1", "Jag2", "Dll1", "Dll3", "Dll4"), group.by = "stage")
dev.off()

png("dot_RPC_target.png")
DotPlot(object = pbmc_RPC, features = c("Hes1", "Hes5", "Hey1", "Hey2"), group.by = "stage")
dev.off()

# sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.1.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.1.3         Seurat_4.9.9.9060       SeuratObject_4.9.9.9091
# [4] sp_2.1-1                dplyr_1.1.4            
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
# [4] rlang_1.1.2            magrittr_2.0.3         RcppAnnoy_0.0.21      
# [7] matrixStats_1.0.0      ggridges_0.5.4         compiler_4.3.1        
# [10] spatstat.geom_3.2-7    png_0.1-8              vctrs_0.6.4           
# [13] reshape2_1.4.4         stringr_1.5.0          pkgconfig_2.0.3       
# [16] fastmap_1.1.1          ellipsis_0.3.2         labeling_0.4.3        
# [19] utf8_1.2.4             promises_1.2.1         purrr_1.0.2           
# [22] jsonlite_1.8.7         goftest_1.2-3          later_1.3.1           
# [25] spatstat.utils_3.0-3   irlba_2.3.5.1          parallel_4.3.1        
# [28] cluster_2.1.4          R6_2.5.1               ica_1.0-3             
# [31] stringi_1.7.12         RColorBrewer_1.1-3     spatstat.data_3.0-1   
# [34] reticulate_1.34.0      parallelly_1.36.0      lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.0.11            tensor_1.5            
# [40] future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.1     
# [43] httpuv_1.6.11          Matrix_1.6-1.1         splines_4.3.1         
# [46] igraph_1.5.1           tidyselect_1.2.0       rstudioapi_0.15.0     
# [49] abind_1.4-5            spatstat.random_3.1-6  codetools_0.2-19      
# [52] miniUI_0.1.1.1         spatstat.explore_3.2-3 listenv_0.9.0         
# [55] lattice_0.21-9         tibble_3.2.1           plyr_1.8.9            
# [58] withr_2.5.2            shiny_1.7.5.1          ROCR_1.0-11           
# [61] Rtsne_0.16             future_1.33.0          fastDummies_1.7.3     
# [64] survival_3.5-7         polyclip_1.10-6        fitdistrplus_1.1-11   
# [67] pillar_1.9.0           KernSmooth_2.23-22     plotly_4.10.3         
# [70] generics_0.1.3         RcppHNSW_0.5.0         ggplot2_3.4.4         
# [73] munsell_0.5.0          scales_1.2.1           globals_0.16.2        
# [76] xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2        
# [79] tools_4.3.1            data.table_1.14.8      RSpectra_0.16-1       
# [82] RANN_2.6.1             leiden_0.4.3           dotCall64_1.1-0       
# [85] cowplot_1.1.1          grid_4.3.1             tidyr_1.3.0           
# [88] colorspace_2.1-0       nlme_3.1-163           cli_3.6.1             
# [91] spatstat.sparse_3.0-2  spam_2.9-1             fansi_1.0.5           
# [94] viridisLite_0.4.2      uwot_0.1.16            gtable_0.3.4          
# [97] digest_0.6.33          progressr_0.14.0       ggrepel_0.9.4         
# [100] farver_2.1.1           htmlwidgets_1.6.2      htmltools_0.5.6.1     
# [103] lifecycle_1.0.4        httr_1.4.7             mime_0.12             
# [106] MASS_7.3-60   