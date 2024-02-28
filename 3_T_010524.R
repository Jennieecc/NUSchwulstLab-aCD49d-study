Tcells <- subset(s, idents=c("T"))
Tcells<-FindVariableFeatures(Tcells)
Tcells <- RunPCA(object = Tcells, verbose = FALSE, assay = "integrated")
Tcells <- RunUMAP(Tcells, dims = 1:30, verbose = FALSE)
Tcells <- FindNeighbors(object = Tcells, dims = 1:30)

#Projecting scRNA-seq data onto a reference map of tumour-infiltrating lymphocytes (TILs)
library(remotes)
library(ProjecTILs)
library(TILPRED)
library(dplyr)
library(tibble)
library(Seurat)
ref <- load.reference.map()
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")

#use ProjecTIL to predict T states
querydata <- Tcells
query.projected <- make.projection(querydata, ref=ref, skip.normalize = TRUE,scGate_model = NULL, query.assay = "RNA")
plot.projection(ref, query.projected)
query.projected <- cellstate.predict(ref=ref, query=query.projected)
plot.statepred.composition(ref, query.projected,metric = "Percent")
plot.states.radar(ref, query=query.projected_clonal, min.cells = 20)
Idents(query.projected)<-query.projected$functional.cluster
markers <- c("Ccr7", "Cd4", "Cd8a", "Pdcd1", "Tox","Tcf7","Lef1","Il7r","Gzmb", "Gzmk","Izumo1r",
           "Ifng", "Foxp3")
query.projected <- RunUMAP(query.projected, dims = 1:30, verbose = FALSE)
query.projected <- FindNeighbors(object = query.projected, dims = 1:30)
query.projected <- FindClusters(object = query.projected, resolution = 1.2) 
Idents(query.projected)<-query.projected$seurat_clusters
DimPlot(query.projected)
DefaultAssay(object = query.projected_clonal) <- "RNA"
query.projected <- NormalizeData(query.projected, verbose = FALSE)
VlnPlot(query.projected, features = markers_Th, stack = T, flip = T,sort = TRUE, same.y.lims = TRUE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#Combine single cell RNAseq with TCRseq
Idents(query.projected)<-query.projected$clonal
DimPlot(query.projected,group.by = "functional.cluster", cols = refCols)

query.projected_clonal <- subset(query.projected, clonal %in% c("NC", "C"))

pt1 <- table(query.projected_clonal$orig.ident, query.projected_clonal$clonal)
pt1 <- as.data.frame(pt1)
pt1$Var1 <- as.character(pt$Var1)
ClonalvsNonclonal1<-ggplot(pt1, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("green", "darkblue"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
Idents(query.projected_clonal)<-query.projected_clonal$functional.cluster
DimPlot(query.projected_clonal, cols = refCols)
Idents(query.projected_clonal)<-query.projected_clonal$clonal
DimPlot(query.projected_clonal, cols = c("darkblue","green"))

Idents(query.projected_clonal)<-query.projected_clonal$functional.cluster
query.projected_clonal_CD8EF<-subset(query.projected_clonal, idents=c("CD8_EffectorMemory"))
pt2 <- table(query.projected_clonal_CD8EF$orig.ident, query.projected_clonal_CD8EF$clonal)
pt2 <- as.data.frame(pt2)

ClonalvsNonclonal_CD8EF<-ggplot(pt2, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("green", "darkblue"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

Idents(query.projected_clonal)<-query.projected_clonal$functional.cluster
query.projected_clonal_CD8NL<-subset(query.projected_clonal, idents=c("CD8_NaiveLike"))
pt3 <- table(query.projected_clonal_CD8NL$orig.ident, query.projected_clonal_CD8NL$clonal)
pt3 <- as.data.frame(pt3)
ClonalvsNonclonal_CD8NL<-ggplot(pt3, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("green", "darkblue"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

ClonalvsNonclonal1+ClonalvsNonclonal_CD8EF+ClonalvsNonclonal_CD8NL

# Divide clonality by the frequency of clonotypes
query.projected_clonaltype<-query.projected_clonal
query.projected_clonaltype$clonal <- "NA"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency == 1] <- "Single"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >1 &  query.projected_clonaltype$frequency <10] <- "Small"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >=35 &  query.projected_clonaltype$frequency <45] <- "Medium"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >=60 & query.projected_clonaltype$frequency <=101] <- "Large"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency > 101] <- "Hyperexpanded"

table(query.projected_clonaltype@meta.data$orig.ident,query.projected_clonaltype@meta.data$frequency)

Idents(query.projected_clonaltype)<-query.projected_clonaltype$clonal

DimPlot(query.projected_clonaltype, label = FALSE ,repel = TRUE,cols = c("grey","dodgerblue3","lightblue","darkgoldenrod1","red"))

pt <- table(query.projected_clonaltype$orig.ident, query.projected_clonaltype$clonal)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 13) +
  geom_col(position = "fill", width = 0.7) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("red","darkgoldenrod1","dodgerblue3","lightblue","lightgrey"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

pt <- table(query.projected_clonaltype$functional.cluster, query.projected_clonaltype$clonal)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("red","darkgoldenrod1","lightblue"))

# scoring of top expressed genes in the clonally expanded CD8 EF cells 
query.projected_C_CD8EF <- AddModuleScore(query.projected_C_CD8EF, features = c("Tox"), name = "Tox")
query.projected_C_CD8EF <- AddModuleScore(query.projected_C_CD8EF, features = c("Rbm3"), name = "Rbm3")
query.projected_C_CD8EF <- AddModuleScore(query.projected_C_CD8EF, features = c("Ifngr1"), name = "Ifngr1")

  VlnPlot(query.projected_C_CD8EF,c("Ifngr1"),group.by = "orig.ident", pt.size = 0.1, adjust = 0)+
  geom_violin()+
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6,col = "black")

#Create df for Rbm3 and Tox
  metadata_ClonalCD8EF <- query.projected_C_CD8EF@meta.data
  mg_enri=metadata_ClonalCD8EF[,c(1,34)]
  rownames(mg_enri)=NULL
  mg_enri$orig.ident=as.factor(mg_enri$orig.ident)
  mg_enri_group=levels(mg_enri$orig.ident)
  str(mg_enri)
  
  mg_enri.pairs <- combn(seq_along(mg_enri_group), 2, simplify=FALSE, FUN=function(i)mg_enri_group[i])
  
  library(ggpubr)
  p1 <- ggplot(mg_enri, aes(x = orig.ident, y = Ifngr11))+geom_violin(trim=FALSE, fill = "palegreen") 
  
  ####Pairwise comparision using non-parametric test (Wilcoxon test)
  p1=p1+stat_compare_means(comparisons = mg_enri.pairs)#, label = "p.signif"


  
  # Find DEGs 
  Idents(query.projected_clonal)<-query.projected_clonal$clonal
  query.projected_C<-subset(query.projected_clonal, idents=c("C"))
  query.projected_NC<-subset(query.projected_clonal, idents=c("NC"))
  
  Idents(query.projected_C)<-query.projected_C$functional.cluster
  Idents(query.projected_NC)<-query.projected_NC$functional.cluster
  
  #CD8EF
  query.projected_C_CD8EF<-subset(query.projected_C, idents=c("CD8_EffectorMemory"))
  query.projected_NC_CD8EF<-subset(query.projected_NC, idents=c("CD8_EffectorMemory"))
  
  Idents(query.projected_C_CD8EF)<-query.projected_C_CD8EF$orig.ident
  Idents(query.projected_NC_CD8EF)<-query.projected_NC_CD8EF$orig.ident
  
  DefaultAssay(object = query.projected_C_CD8EF) <- "RNA"
  query.projected_C_CD8EF <- NormalizeData(query.projected_C_CD8EF, verbose = FALSE)
  
  DefaultAssay(object = query.projected_NC_CD8EF) <- "RNA"
  query.projected_NC_CD8EF <- NormalizeData(query.projected_NC_CD8EF, verbose = FALSE)
  
  
  degs <-FindMarkers(object = query.projected_C_CD8EF,
                      ident.1 = "ATBI_aCD49d",
                      ident.2 = "ATBI_Isotype",
                      test.use = "MAST",
                      logfc.threshold = -Inf,
                      min.pct = 0.5,
                      assay = "RNA"
  ) 
  
  degs2 <-FindMarkers(object = query.projected_NC_CD8NL,
                     ident.1 = "ATBI_aCD49d",
                     ident.2 = "ATBI_Isotype",
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.5,
                     assay = "RNA"
  )  
  
  # Remove ribosomal, mitochondrial, and HLA genes
  degs <- degs[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(degs)),]
  degs <- degs[!grepl("AU", rownames(degs)),]
  degs <- degs[!grepl("24", rownames(degs)),]
  
  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  degs$genes <- rownames(degs)
  
  # Write out results
  write.csv(degs2, "CD8_NC_EF_ATBIaCD49d_vs_ATBI_Isotype_degs.csv")
  #repeat for all the T cell subsets
  deg_merge <- read.csv("~_CD8_T_Merge_ATBIaCD49d_vs_ATBI_Isotype_degs.csv")


  # Create volcano plot
  options(ggrepel.max.overlaps = 10)
#use 00_volcanoMerge_function 
  volcano_plot(deg_merge, title = "CD8+ T Cells with aCD49d Ab 
  Treatment in Aged Mice PTBI",
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  
 
# Run DE on all celltypes
#cell_types <- c("CD8_EffectorMemory", "CD8_NaiveLike", "Th1","CD4_NaiveLike")
#mclapply(cell_types, run_de, s = query.projected_aged_clonal, mc.cores = 3)
#mclapply(cell_types, run_de, s = query.projected_aged_nonclonal, mc.cores = 3)


# Initialize sig gene lists
sig_genes_clonal_ls <- list()

# Create list with sig genes for each cell type
cell_types <- c("CD8EA_C","CD8EA_NC","CD8EF_C","CD8EF_NC","CD8NL_C","CD8NL_NC","Th1_C","Th1_NC")
  
for (cell_type in cell_types) {
  # Define cell type label
  cell_type_label <- gsub("/", "", cell_types)
  cell_type_label <- gsub(" ", "_", cell_type_label)

  # Load in degs
 #output_dir <- "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/TCRSeuratmerge/"
  degs_c <- read.csv(paste0(output_dir, cell_type_label, "_ATBIaCD49d_vs_ATBI_Isotype_degs.csv"))
  
  # Identify sig genes
#  padj.thresh <- 0.01
 # lfc.thresh <- 0.25
  sig_genes_c <- degs_c[which(degs_c$BH <= padj.thresh & abs(degs_c$avg_log2FC) >= lfc.thresh),]
  
  # Add sig genes to list
  sig_genes_clonal_ls[[cell_type]] <- sig_genes_c$X
}

# Find number of genes in each set
num_degs_c <- sapply(sig_genes_clonal_ls, length)
sig_genes_clonal_ls <- sig_genes_clonal_ls[order(-num_degs_c)]
sig_genes_clonal_ls <- sig_genes_clonal_ls[ lapply(sig_genes_clonal_ls, length) > 0 ]

# Define colors
color_vector_c <- mapvalues(
  names(sig_genes_clonal_ls),
  from = c("CD8EA_C","CD8EA_NC","CD8EF_C","CD8EF_NC","CD8NL_C","CD8NL_NC","Th1_C","Th1_NC"),
  to = c("#A58AFF", "#A58AFF", "#53B400","#53B400", "#edbe2a","#edbe2a","palegreen","palegreen")
)

# Create upset plot and export
pdf(file = paste0(output_dir, "/clonal_CI_vs_HC_upset_celltype.pdf"))
upset(
  fromList(sig_genes_clonal_ls),
  nsets = length(sig_genes_clonal_ls),
  order.by = "freq",
  sets.bar.color = color_vector_c,
  text.scale = 2, nintersects = 15
)
dev.off()



