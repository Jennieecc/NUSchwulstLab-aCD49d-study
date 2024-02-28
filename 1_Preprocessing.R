####################################
# data preprocessing and analysis  #
####################################

# Load packages
library(parallel)
library(future)
library(Seurat)
library(Matrix)
options(future.globals.maxSize=256*1024^3)
plan(multiprocess)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(tidyverse)
library(purrr)

# Analysis parameters
anchor_dims <- 30            # number of anchor dimensions used for biological replicates integration
pca_dims <- 30               # number of PCA dimensions to compute and use in tSNE, UMAP                             
umap_n_neighbors <- 30       # UMAP parameter,                                                                    
clustering_resolution <- 0.8 # Resolution parameter for Seurat clustering
n_features <- 2000

#SoupX to remove cell free mRNA contamination 
library(SoupX)
sc1 = load10X('/fsmresfiles.fsm.northwestern.edu/fsmresfiles/Surgery/Schwulst_Lab/Jennie_Chen/Single_Cell_12122022/S9_H/outs')
sc2 = load10X('/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/2/outs')
sc3 = load10X('/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/3/outs')
sc4 = load10X('/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/4/outs')

sc1 = autoEstCont(sc1)
sc2 = autoEstCont(sc2)
sc3 = autoEstCont(sc3)
sc4 = autoEstCont(sc4)

YTBI_Isotype = adjustCounts(sc1)
YTBI_aCD49d = adjustCounts(sc2)
ATBI_Isotype = adjustCounts(sc3)
ATBI_aCD49d = adjustCounts(sc4)

YTBI_Isotype = CreateSeuratObject(YTBI_Isotype)
YTBI_aCD49d = CreateSeuratObject(YTBI_aCD49d)
ATBI_Isotype = CreateSeuratObject(ATBI_Isotype)
ATBI_aCD49d = CreateSeuratObject(ATBI_aCD49d)

YTBI_Isotype$orig.ident<-"YTBI_Isotype"
YTBI_aCD49d$orig.ident<-"YTBI_aCD49d"
ATBI_Isotype$orig.ident <- "ATBI_Isotype"
ATBI_aCD49d$orig.ident <- "ATBI_aCD49d"

# Analyze percentage of mitochondrial genes in cells
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
YTBI_Isotype[["percent.mt"]] <- PercentageFeatureSet(YTBI_Isotype, pattern = "^mt-")
YTBI_aCD49d[["percent.mt"]] <- PercentageFeatureSet(YTBI_aCD49d, pattern = "^mt-")
ATBI_Isotype[["percent.mt"]] <- PercentageFeatureSet(ATBI_Isotype, pattern = "^mt-")
ATBI_aCD49d[["percent.mt"]] <- PercentageFeatureSet(ATBI_aCD49d, pattern = "^mt-")

#metrics to measure are cell counts, UMI counts, genes detected/cell, UMIs vs genes detected, Mito ratio, and novelty
YTBI_Isotype$log10GenesPerUMI <- log10(YTBI_Isotype$nFeature_RNA)/log10(YTBI_Isotype$nCount_RNA)
YTBI_aCD49d$log10GenesPerUMI <- log10(YTBI_aCD49d$nFeature_RNA)/log10(YTBI_aCD49d$nCount_RNA)
ATBI_Isotype$log10GenesPerUMI <- log10(ATBI_Isotype$nFeature_RNA)/log10(ATBI_Isotype$nCount_RNA)
ATBI_aCD49d$log10GenesPerUMI <- log10(ATBI_aCD49d$nFeature_RNA)/log10(ATBI_aCD49d$nCount_RNA)

# Create metadata dataframe
metadata_ATBI_aCD49d <- ATBI_aCD49d@meta.data

# Add cell IDs to metadata
metadata_ATBI_aCD49d$cells <- rownames(metadata_ATBI_aCD49d)

# Rename columns
metadata_ATBI_aCD49d <- metadata_ATBI_aCD49d %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
ATBI_aCD49d@meta.data <- metadata_ATBI_aCD49d

#visualize/assess quality metrics
# Visualize the number of cell counts per sample
ATBI_aCD49d@meta.data%>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar() +
  theme_classic() +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("ATBI_aCD49d NCells")

# Visualize the number UMIs/transcripts per cell
ATBI_aCD49d@meta.data %>% 
  ggplot(aes(x=nUMI)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1400)

# ATBI_Isotype the distribution of genes detected per cell via histogram
ATBI_aCD49d@meta.data %>% 
  ggplot(aes( x=nGene)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 660)

# Visualize the distribution of genes detected per cell via boxplot
ATBI_aCD49d@meta.data %>% 
  ggplot(aes(y=log10(nGene))) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
ATBI_aCD49d@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1400) +
  geom_hline(yintercept = 660) 

# Visualize the distribution of mitochondrial gene expression detected per cell
ATBI_aCD49d@meta.data %>% 
  ggplot(aes(x=percent.mt)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
ATBI_aCD49d@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Filter out low quality cells using selected thresholds - these will change with experiment
ATBI_aCD49d_filtered <- subset(x = ATBI_aCD49d, 
                               subset= (nUMI >= 1400) & 
                                 (nGene >= 660) & 
                                 (log10GenesPerUMI > 0.80) & 
                                 (percent.mt < 10))




#Gene-level filtering
# Extract counts
counts <- GetAssayData(object = YTBI_Isotype_filtered, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
YTBI_Isotype_filtered <- CreateSeuratObject(filtered_counts, meta.data = YTBI_Isotype_filtered@meta.data)


# Extract counts
counts <- GetAssayData(object = YTBI_Isotype_filtered, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
YTBI_Isotype_filtered <- CreateSeuratObject(filtered_counts, meta.data = YTBI_Isotype_filtered@meta.data)

#explore sources of unwanted variation
# Normalize the counts
seurat_phase_aged <- NormalizeData(ATBI_aCD49d_filtered)
cell_cycle_genes <- readLines(con = "~/Desktop/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
s.genes <- cell_cycle_genes[1:43]
g2m.genes <- cell_cycle_genes[44:97]
seurat_phase_aged <- CellCycleScoring(seurat_phase_aged, 
                                      g2m.features = g2m.genes, 
                                      s.features = s.genes)
# Identify the most variable genes
seurat_phase_aged <- FindVariableFeatures(seurat_phase_aged, 
                                          selection.method = "vst",
                                          nfeatures = 2000, 
                                          verbose = FALSE)

# Scale the counts
seurat_phase_aged <- ScaleData(seurat_phase_aged)

# Perform PCA
seurat_phase_aged <- RunPCA(seurat_phase_aged)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase_aged,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

#We do not see large differences due to cell cycle phase in both young merge and old merge. Therefore, we would not regress out the variation due to cell cycle

#within group variation
all_obj <- merge(x=YTBI_Isotype_filtered,y=c(YTBI_aCD49d_filtered,ATBI_Isotype_filtered, ATBI_aCD49d_filtered), add.cell.ids = c("YTBI_Isotype","YTBI_aCD49d", "ATBI_Isotype", "ATBI_aCD49d"))

all_obj$orig.ident
all_obj<-NormalizeData(all_obj)
mean_matrix <- AverageExpression(all_obj, slot = "data", group.by = "orig.ident")$"RNA"
cor_matrix <- cor(mean_matrix) %>% reshape2::melt()

ggplot(cor_matrix, aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+viridis::scale_fill_viridis()+ theme_classic(base_size = 24)

#Integration as young samples and aged samples are from different batch
#using SCTransform objects as input run on Quest (R)
YTBI_Isotype_SCT <- SCTransform(YTBI_Isotype_filtered,ncells = 5000, conserve.memory = TRUE)
YTBI_aCD49d_SCT <- SCTransform(YTBI_aCD49d_filtered,ncells = 5000, conserve.memory = TRUE)
ATBI_Isotype_SCT <- SCTransform(ATBI_Isotype_filtered,ncells = 5000, conserve.memory = TRUE)
ATBI_aCD49d_SCT <- SCTransform(ATBI_aCD49d_filtered,ncells = 5000, conserve.memory = TRUE)


saveRDS(YTBI_Isotype_SCT, file = "YTBI_Isotype_SCT.rds")
saveRDS(YTBI_aCD49d_SCT, file = "YTBI_aCD49d_SCT.rds")
saveRDS(ATBI_Isotype_SCT, file = "ATBI_Isotype_SCT.rds")
saveRDS(ATBI_aCD49d_SCT, file = "ATBI_aCD49d_filtered.rds")

ATBI_aCD49d_SCT <- readRDS(file = "/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/ATBI_aCD49d_SCT.rds")
YTBI_aCD49d_SCT <- readRDS(file = "/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/YTBI_aCD49d_SCT.rds")
ATBI_Isotype_SCT <- readRDS(file = "/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/ATBI_Isotype_SCT.rds")
YTBI_Isotype_SCT <- readRDS(file = "/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/YTBI_Isotype_SCT.rds")

#samples_integrated_2wk_8hr_M.rds <- readRDS(file = "/Users/zhangyingchen/Desktop/Single_Cell_12122022/samples_integrated_2wk_8hr_M.rds")

#aged_merge_SCF <- merge(x=split_filtered_aged$sham,y=split_filtered_aged$TBI, add.cell.ids = c("aged_sham","aged_TBI"))
#aged_merge_SCF <- SCTransform(aged_merge_SCF, vars.to.regress = c("mitoRatio"))

#var_genes_young <- Reduce(intersect, list(VariableFeatures(split_filtered_young$sham), VariableFeatures(split_filtered_young$TBI)))
#var_genes_aged <- intersect(VariableFeatures(split_filtered_aged$sham),VariableFeatures(split_filtered_aged$TBI))
#VariableFeatures(young_merge_SCF) <- var_genes_young
#VariableFeatures(aged_merge_SCF) <- var_genes_aged

features <- SelectIntegrationFeatures(object.list = c(YTBI_Isotype_SCT, YTBI_aCD49d_SCT, ATBI_Isotype_SCT,ATBI_aCD49d_SCT), nfeatures = 3000)
merge.list <- PrepSCTIntegration(object.list = c(YTBI_Isotype_SCT, YTBI_aCD49d_SCT, ATBI_Isotype_SCT,ATBI_aCD49d_SCT), anchor.features = 3000)

library(future)
plan("multiprocess", workers = 4)
samples_anchors <- FindIntegrationAnchors(object.list = merge.list, 
                                          dims = 1:30,
                                          anchor.features = features,
                                          normalization.method="SCT")

samples_integrated_CD49d <- IntegrateData(anchorset = samples_anchors, dims = 1:30, normalization.method = "SCT")
samples_integrated_CD49d<- FindVariableFeatures(samples_integrated_CD49d)
samples_integrated_CD49d <- RunPCA(object = samples_integrated_CD49d, verbose = FALSE)
samples_integrated_CD49d <- RunUMAP(samples_integrated_CD49d, dims = 1:30, verbose = FALSE)
samples_integrated_CD49d <- FindNeighbors(object = samples_integrated_CD49d, dims = 1:30)
samples_integrated_CD49d <- FindClusters(object = samples_integrated_CD49d, resolution = 0.8) 
#samples_integrated_CD49d$sample <- paste(samples_integrated_CD49d$orig.ident)#, all_obj$sample, sep = "_"
DimPlot(samples_integrated_CD49d,label = TRUE, split.by = "orig.ident", ncol = 2)
DimPlot(samples_integrated_CD49d,label = TRUE)
saveRDS(samples_integrated_CD49d, file = "/Users/zhangyingchen/Desktop/AvY antiCD49d study/scRNA/samples_integrated_CD49d.rds")

pt <- table(samples_integrated_CD49d$pred.celltypes, samples_integrated_CD49d$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())
  scale_fill_manual(values=refCols)

#SingleR annotation
library(SingleR)
library(celldex)
moref <- MouseRNAseqData()
pred.samples_integrated <- SingleR(test = samples_integrated_CD49d@assays$RNA@data, ref = moref, labels = moref$label.main)
table(pred.samples_integrated$labels)
samples_integrated_CD49d$pred.celltypes <- pred.samples_integrated$labels

Idents(samples_integrated_CD49d)=samples_integrated_CD49d$pred.celltypes
DimPlot(s, reduction="umap", label = TRUE)
DimPlot(samples_integrated_CD49d, split.by = c("orig.ident"), ncol = 2,label = TRUE)
s <- subset(samples_integrated_CD49d, idents=c("T cells","B cells","Granulocytes","Macrophages",
                                               "Microglia","Monocytes","NK cells"))

s <- RunPCA(object = s, verbose = FALSE)
s <- RunUMAP(s, dims = 1:30, verbose = FALSE)
s <- FindNeighbors(object = s, dims = 1:30)
s <- FindClusters(object = s, resolution = 0.8) 

#Main figures
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
FeaturePlot(object = s,features =  c("P2ry12","Tmem119","Hexb","C1qa","C1qb","Ctss"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
T_cell <- c("Cd8a","Il7r","Cd3d","Cd4")
FeaturePlot(object = s, features = c("Cd8a","Cd4","Cd3d","Foxp3"), order = TRUE, reduction = "umap", min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
B_cell <- c("Cd79a","Cd79b","Cd74","Ly6d","Ms4a1")
FeaturePlot(object = s, features = B_cell, order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
NK_cell <- c("Ncr1","Klre1") 
FeaturePlot(object = s, features = c("Ncr1","Klre1"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
#CAM in cluster 6&25
CAM <- c("Lyve1","Cd163","Siglec1")
FeaturePlot(object = s, features = c("Lyve1","Cd163","Siglec1","Pf4","Mrc1"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
VlnPlot(samples_integrated,features = c("Lyve1","Siglec1","Cd163"))
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
#CAM in cluster 6&25
#classical monocytes 17,19,20,31
Mo_MΦ <- c("Ccr2","Ly6c2","Tgfbi","Ifitm2","Ifitm3","F13a1","S100a11","S100a6")
FeaturePlot(object = s, features = Mo_MΦ,  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)

FeaturePlot(object = samples_integrated_8hr_MvsF, features = c("Ccr2","Ly6c2","Ly6c1"),  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)

#UK cluster
UK <- c("Meg3","Clu")
FeaturePlot(object = samples_integrated_2wk_8hr_M.rds, features = UK,  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
#CAM in cluster 6&25
#classical monocytes 17,19,20,31
#neutrophil cluster 7&13&26
NP <- c("Cxcr2", "Csf3r","Ly6g","S100a9","S100a8")
FeaturePlot(object = s, features = c("Cxcr2", "Csf3r","Ly6g","S100a9","S100a8"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
#CAM in cluster 6&25
#classical monocytes 17,19,20,31
#neutrophil cluster 7&13&26
#Nonclassical monocytes cluster 21
NonMo <- c("Spn")
FeaturePlot(object = s, features = c("Spn", "Ace"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
VlnPlot(samples_integrated,features = c("Spn", "Ace"))

#Dentritic cells 
FeaturePlot(object = s, features = c("Flt3","Zbtb46","Clec9a"),  order = TRUE, reduction = "umap", ncol=1, min.cutoff = 'q10', label = TRUE)

DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)
all.markers.RNA <- FindAllMarkers(s,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

background <- FindAllMarkers(samples_integrated_2wk_8hr_M.rds,assay="RNA",slot="data", min.pct = 0.25, logfc.threshold = 0)

top15_comb <- all.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)


#DE analysis
DefaultAssay(samples_integrated) <- "RNA"
samples_integrated <- NormalizeData(samples_integrated, verbose = FALSE)
all.markers.RNA <- FindAllMarkers(samples_integrated,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- all.markers.RNA %>% group_by(cluster) %>% top_n(n = 10, avg_log2FC)

allcells_referece <-FindAllMarkers(samples_integrated_CD49d, only.pos = TRUE, min.pct =  0.25, logfc.threshold = 0)

write.csv(allcells_referece, 
          file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/allcells_referece.csv", 
          quote = FALSE, 
          row.names = FALSE)

top15_comb <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)
# Save markers to file
write.csv(top15_comb, 
          file = "/Users/zhangyingchen/Desktop/top15_comb.csv", 
          quote = FALSE, 
          row.names = FALSE)

# Figure 
gene_panel <- c("Cd8a","Il7r","Cd3d","Cd4","P2ry12","Tmem119","C1qa","Ctss","Mertk","Cd79a",
                "Ly6d","Ms4a1","Ttr","Lyve1","Cd163","Siglec1","Krt18","Slc2a1","Ncr1","Klre1","Ly6c1","Ly6c2", "Ccr2","Cxcr2","Ly6g","S100a8", 
                "Ifitm2","Ifitm3","S100a6","Tgfbi")


DotPlot(samples_integrated, features = gene_panel) + 
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "dark blue",high= "red")

VlnPlot(samples_integrated, features = gene_panel, flip = T, stack = T) 

#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
#T cell cluster 0,16,23,24,27,28
#B cell cluster 3&22
#NK cluster cluster 14
#CAM in cluster 6&25
#classical monocytes 17,19,20,31
#neutrophil cluster 7&13&26
#Nonclassical monocytes cluster 21

s$cell_type_8_groups <- plyr::mapvalues(Idents(s), from=c(0:31), 
                                        to=c("T", "MG", "MG", "B", "MG", "MG","MΦ","NP",
                                             "MG","MG","MG","MG","MG","NP","NK","MG","T",
                                             "Mo","MG","Mo","Mo","ncMo","B","T","T","MΦ","NP",
                                             "T","T","B","NP","Mo"))

Idents(s)=s$cell_type_8_groups

main_col <- c("#B0C4DE","#F08080","#9ACD32", "#FFE4B5","#FABF00","#6495ED","#C2B4FC","#DFA5F2")
DimPlot(s, label = FALSE,cols = main_col,pt.size = 0.8)
DimPlot(s, split.by = "orig.ident",label = TRUE, ncol = 2,cols = main_col)
gene_panel <- c("Cd8a","Il7r","Cd3d","Cd4","P2ry12","Tmem119","Mertk","Sall1","Sparc","Cd79a",
                "Ly6d","Ms4a1","Lyve1","Cd163","Siglec1","Cxcr2","Ly6g","S100a8","Ncr1","Klre1","Ly6c2", "Ccr2",
                "Ifitm2","Ifitm3","S100a6","Tgfbi","Ace","Spn")
#"Ttr", "Krt18","Slc2a1", what are these?
DotPlot(s, features = gene_panel) + 
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "dark blue",high= "red")

pt <- table(Idents(s), s$orig.ident)
pt <- as.data.frame(pt)
pt_ATBI_aCD49d = pt[c(1:8),]
pt_ATBI_Isotype = pt[c(9:16),]
pt_YTBI_aCD49d = pt[c(17:24),]
pt_YTBI_Isotype = pt[c(25:32),]

pt_ATBI_aCD49d$Var1 <- as.character(pt_ATBI_aCD49d$Var1)
pt_ATBI_Isotype$Var1 <- as.character(pt_ATBI_Isotype$Var1)
pt_YTBI_aCD49d$Var1 <- as.character(pt_YTBI_aCD49d$Var1)
pt_YTBI_Isotype$Var1 <- as.character(pt_YTBI_Isotype$Var1)

library(RColorBrewer)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values = brewer.pal(10, "Paired"))

saveRDS(s, file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/s.rds")

#cellchat
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

cellchat <- createCellChat(object = s,
                           group.by = "cell_type_8_groups")

# Define ligand-receptor pair database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

# Preprocessing
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculate communication probability and network
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

#calculate communication probability and netwirj
# Define celltype colors

# Visualize outgoing vs incoming signal
netAnalysis_signalingRole_scatter(cellchat,color.use = main_col,dot.size = c(2, 12),
                                  font.size = 14, do.label = T, label.size = 5,
                                  title = "Cell-Cell Interaction Strength",
                                  font.size.title = 15) 



#Seperate cellchat
# Subset for disease
s_ATBI_iso <- subset(s, orig.ident == "ATBI_Isotype")
s_YTBI_iso <- subset(s, orig.ident == "YTBI_Isotype")
s_ATBI_aCD49d <- subset(s, orig.ident == "ATBI_aCD49d")
s_YTBI_aCD49d <- subset(s, orig.ident == "YTBI_aCD49d")

# Run cellchat
cellchat_ATBI_iso <- createCellChat(object = s_ATBI_iso,
                              group.by = "cell_type_8_groups")
cellchat_YTBI_iso <- createCellChat(object = s_YTBI_iso,
                              group.by = "cell_type_8_groups")
cellchat_ATBI_aCD49d <- createCellChat(object = s_ATBI_aCD49d,
                                    group.by = "cell_type_8_groups")
cellchat_YTBI_aCD49d <- createCellChat(object = s_YTBI_aCD49d,
                                    group.by = "cell_type_8_groups")

# Set the used database in the object
cellchat_ATBI_iso@DB <- CellChatDB
cellchat_YTBI_iso@DB <- CellChatDB
cellchat_ATBI_aCD49d@DB <- CellChatDB
cellchat_YTBI_aCD49d@DB <- CellChatDB


# Preprocessing
cellchat_ATBI_iso <- subsetData(cellchat_ATBI_iso)
cellchat_YTBI_iso <- subsetData(cellchat_YTBI_iso)
cellchat_ATBI_iso <- identifyOverExpressedGenes(cellchat_ATBI_iso)
cellchat_ATBI_iso <- identifyOverExpressedInteractions(cellchat_ATBI_iso)
cellchat_YTBI_iso <- identifyOverExpressedGenes(cellchat_YTBI_iso)
cellchat_YTBI_iso <- identifyOverExpressedInteractions(cellchat_YTBI_iso)

cellchat_ATBI_aCD49d <- subsetData(cellchat_ATBI_aCD49d)
cellchat_YTBI_aCD49d <- subsetData(cellchat_YTBI_aCD49d)
cellchat_ATBI_aCD49d <- identifyOverExpressedGenes(cellchat_ATBI_aCD49d)
cellchat_ATBI_aCD49d <- identifyOverExpressedInteractions(cellchat_ATBI_aCD49d)
cellchat_YTBI_aCD49d <- identifyOverExpressedGenes(cellchat_YTBI_aCD49d)
cellchat_YTBI_aCD49d <- identifyOverExpressedInteractions(cellchat_YTBI_aCD49d)

# Calculate communication probability and network
cellchat_ATBI_iso <- computeCommunProb(cellchat_ATBI_iso)
cellchat_ATBI_iso <- computeCommunProbPathway(cellchat_ATBI_iso)
cellchat_ATBI_iso <- aggregateNet(cellchat_ATBI_iso)
cellchat_ATBI_iso <- netAnalysis_computeCentrality(cellchat_ATBI_iso, slot.name = "netP")
cellchat_YTBI_iso <- computeCommunProb(cellchat_YTBI_iso)
cellchat_YTBI_iso <- computeCommunProbPathway(cellchat_YTBI_iso)
cellchat_YTBI_iso <- aggregateNet(cellchat_YTBI_iso)
cellchat_YTBI_iso <- netAnalysis_computeCentrality(cellchat_YTBI_iso, slot.name = "netP")

cellchat_ATBI_aCD49d <- computeCommunProb(cellchat_ATBI_aCD49d)
cellchat_ATBI_aCD49d <- computeCommunProbPathway(cellchat_ATBI_aCD49d)
cellchat_ATBI_aCD49d <- aggregateNet(cellchat_ATBI_aCD49d)
cellchat_ATBI_aCD49d <- netAnalysis_computeCentrality(cellchat_ATBI_aCD49d, slot.name = "netP")
cellchat_YTBI_aCD49d <- computeCommunProb(cellchat_YTBI_aCD49d)
cellchat_YTBI_aCD49d <- computeCommunProbPathway(cellchat_YTBI_aCD49d)
cellchat_YTBI_aCD49d <- aggregateNet(cellchat_YTBI_aCD49d)
cellchat_YTBI_aCD49d <- netAnalysis_computeCentrality(cellchat_YTBI_aCD49d, slot.name = "netP")

netAnalysis_signalingRole_scatter(cellchat_ATBI_iso,color.use = main_col,dot.size = c(2, 12),
                                  font.size = 14, do.label = F, label.size = 4,
                                  title = "Cell-Cell Interaction Strength 
in Aged Isotype Mouse Brains",
                                  font.size.title = 14) 

netAnalysis_signalingRole_scatter(cellchat_YTBI_iso,color.use = main_col,dot.size = c(2, 12),
                                  font.size = 14, do.label = F, label.size = 5,
                                  title = "Cell-Cell Interaction Strength
in Young Isotype Mouse Brains",
                                  font.size.title = 15) 

netAnalysis_signalingRole_scatter(cellchat_ATBI_aCD49d,color.use = main_col,dot.size = c(2, 12),
                                  font.size = 14, do.label = F, label.size = 4,
                                  title = "Cell-Cell Interaction Strength 
in Aged aCD49d Mouse Brains",
                                  font.size.title = 14) 

netAnalysis_signalingRole_scatter(cellchat_YTBI_aCD49d,color.use = main_col,dot.size = c(2, 12),
                                  font.size = 14, do.label = F, label.size = 5,
                                  title = "Cell-Cell Interaction Strength
in Young aCD49d Mouse Brains",
                                  font.size.title = 15) 


# Create cellchat object
netVisual_bubble(cellchat_ATBI_iso,
                 sources.use = 2,
                 targets.use = 1,
                 remove.isolate = FALSE,
                 title.name = "Signaling in Aged Isotype 
  Mouse Brains",
                 font.size =13,
                 font.size.title = 13)
                 #pairLR.use = pairLR.use)
netVisual_bubble(cellchat_YTBI_iso,
                 sources.use = 2,
                 targets.use = 1,
                 title.name = "Signaling in Young Isotype 
  Mouse Brains",  font.size =13,
                 font.size.title = 13,
                 remove.isolate = FALSE)
                 #pairLR.use = pairLR.use)

netVisual_bubble(cellchat_ATBI_aCD49d,
                 sources.use = 2,
                 targets.use = 1,
                 remove.isolate = FALSE,
                 title.name = "Signaling in Aged aCD49d Mouse Brains")

netVisual_bubble(cellchat_YTBI_aCD49d,
                 sources.use = 2,
                 targets.use = 1,
                 title.name = "Signaling in Young aCD49d Mouse Brains",
                 remove.isolate = FALSE)

#merge inididual dataset for changes in signaling
object.list_youngCd49dVSyoungIso <- list(YTBI_Isotype = cellchat_YTBI_iso,  YTBI_aCD49d= cellchat_YTBI_aCD49d)
cellchat3 <- mergeCellChat(object.list_youngCd49dVSyoungIso, add.names = names(object.list_youngCd49dVSyoungIso))


object.list_ageVSyoungIso <- list(YTBI_Isotype = cellchat_YTBI_iso,  ATBI_Isotype= cellchat_ATBI_iso)
cellchat <- mergeCellChat(object.list_ageVSyoungIso, add.names = names(object.list_ageVSyoungIso))

object.list_ageCd49dVSageIso <- list(ATBI_Isotype = cellchat_ATBI_iso, ATBI_aCD49d= cellchat_ATBI_aCD49d)
cellchat_2 <- mergeCellChat(object.list_ageCd49dVSageIso, add.names = names(object.list_ageCd49dVSageIso))

#compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat3, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat3, show.legend = F, group = c(1,2), measure = "weight")
gg3 <- compareInteractions(cellchat_2, show.legend = F, group = c(1,2))
gg4 <- compareInteractions(cellchat_2, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2 +gg3 +gg4

#compare the number of interactions and interactions strength among different cell populations
par(mfrow = c(1,2,3,4), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
netVisual_diffInteraction(cellchat_2, weight.scale = T)
netVisual_diffInteraction(cellchat_2, weight.scale = T, measure = "weight")

#signal changes of MG and T in each comparison
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MG")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_2, idents.use = "MG")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat_2, idents.use = "T")
patchwork::wrap_plots(plots = list(gg3,gg4))

saveRDS(cellchat_ATBI_aCD49d, file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/cellchat_ATBI_aCD49d.rds")
saveRDS(cellchat_ATBI_iso, file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/ccellchat_ATBI_iso.rds")
saveRDS(cellchat_YTBI_aCD49d, file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/cellchat_YTBI_aCD49d.rds")
saveRDS(cellchat_YTBI_iso, file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/cellchat_YTBI_iso.rds")

#Identify dysfunctional signaling by comparing the communication probabilities
netVisual_bubble(cellchat_2, sources.use = c(2,4), targets.use = 1,  remove.isolate = T,
                 font.size = 15,comparison = c(1:2), n.colors = 10,
                 angle.x = 45,thresh = 0.01)

netVisual_bubble(cellchat_2, sources.use =1, targets.use = c(2,4),  remove.isolate = T,
                 font.size = 15,comparison = c(1:2), n.colors = 10,
                 angle.x = 45,thresh = 0.01)
pairLR.use <- extractEnrichedLR(cellchat_2, signaling = c("CCL","TNF","CXCL","MHC-I"))

netVisual_bubble(cellchat_iso, sources.use = 2, targets.use = 1,
                 comparison =  c(1,2), max.dataset = 2,
                 title.name = "Increased signaling with Aging",
                 angle.x = 90, remove.isolate = T,pairLR.use = pairLR.use)

#Comparing DE analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "ATBI_aCD49d"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS and cellchat@var.features$LS.info. 
cellchat_2 <- identifyOverExpressedGenes(cellchat_2, group.dataset = NULL, pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 1) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_2, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_2, net = net, datasets = "ATBI_aCD49d",ligand.logFC = 0.1, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_2, net = net, datasets = "ATBI_Isotype",ligand.logFC = -0.05, receptor.logFC = -0.0)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_2)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_2)

object.list <- list(ATBI_Isotype = cellchat_ATBI_iso, ATBI_aCD49d = cellchat_ATBI_aCD49d)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_2, pairLR.use = pairLR.use.up, sources.use = c(2,4), targets.use = 1, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_2, pairLR.use = pairLR.use.down, sources.use = c(2,4), targets.use = 1, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


pt <- table(cellchat_2@LR$ATBI_aCD49d$LRsig$pathway_name,cellchat_2@LR$ATBI_aCD49d$LRsig$ligand,cellchat_2@LR$ATBI_aCD49d$LRsig$receptor)
pt <- as.data.frame(pt)
