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
sc1 = load10X('~scRNA/1/outs')
sc2 = load10X('~scRNA/2/outs')
sc3 = load10X('~scRNA/3/outs')
sc4 = load10X('~~scRNA/4/outs')

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
saveRDS(samples_integrated_CD49d, file = "~samples_integrated_CD49d.rds")

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
FeaturePlot(object = s,features =  c("P2ry12","Tmem119","Hexb","C1qa","C1qb","Ctss"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#Microglia cluster 1,2,4,5,8,9,10,11,12,15,18
T_cell <- c("Cd8a","Il7r","Cd3d","Cd4")
FeaturePlot(object = s, features = c("Cd8a","Cd4","Cd3d","Foxp3"), order = TRUE, reduction = "umap", min.cutoff = 'q10', label = TRUE)
#T cell cluster 0,16,23,24,27,28
B_cell <- c("Cd79a","Cd79b","Cd74","Ly6d","Ms4a1")
FeaturePlot(object = s, features = B_cell, order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#B cell cluster 3&22
NK_cell <- c("Ncr1","Klre1") 
FeaturePlot(object = s, features = c("Ncr1","Klre1"), order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#NK cluster cluster 14
CAM <- c("Lyve1","Cd163","Siglec1")
FeaturePlot(object = s, features = c("Lyve1","Cd163","Siglec1","Pf4","Mrc1"), order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#CAM in cluster 6&25
Mo_MΦ <- c("Ccr2","Ly6c2","Tgfbi","Ifitm2","Ifitm3","F13a1","S100a11","S100a6")
FeaturePlot(object = s, features = Mo_MΦ,  order = TRUE, reduction = "umap", ncol=3, min.cutoff = 'q10', label = TRUE)
#classical monocytes 17,19,20,31
NP <- c("Cxcr2", "Csf3r","Ly6g","S100a9","S100a8")
FeaturePlot(object = s, features = NP, order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#neutrophil cluster 7&13&26
NonMo <- c("Spn","Ace")
FeaturePlot(object = s, features = NonMo, order = TRUE, reduction = "umap", ncol=2, min.cutoff = 'q10', label = TRUE)
#Nonclassical monocytes cluster 21
#Dentritic cells 
FeaturePlot(object = s, features = c("Flt3","Zbtb46","Clec9a"),  order = TRUE, reduction = "umap", ncol=1, min.cutoff = 'q10', label = TRUE)

DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)
all.markers.RNA <- FindAllMarkers(s,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
background <- FindAllMarkers(s,assay="RNA",slot="data", min.pct = 0.25, logfc.threshold = 0)
top15_comb <- all.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)

# Figure 
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
DotPlot(s, features = gene_panel) + 
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "dark blue",high= "red")
