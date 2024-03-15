## Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("scales")
  library("doMC")
  library("UpSetR")
  library("colorspace")
})
# Set core number for parallel model fitting
registerDoMC(cores = 3)
# Specify thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

Mac <- subset(s, idents=c("MΦ"))
Mac <- RunPCA(object = Mac, verbose = FALSE)
Mac <- RunUMAP(Mac, dims = 1:30, verbose = FALSE)
Mac <- FindNeighbors(object = Mac, dims = 1:30)
Mac <- FindClusters(object = Mac, resolution = 0.8) 
DimPlot(Mac)

Mac.markers.RNA <- FindAllMarkers(Mac,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top15_comb <- Mac.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)

Mac$level <- plyr::mapvalues(Idents(Mac), from=c("0","1","2","3","4","5","6"),
                              to=c("MHCIIhi BAM","MDM","MHCIIlo BAM_1","MHCIIlo BAM_2","MDM","MHCIIlo BAM_2","Ependymal cells"))

Idents(Mac)=Mac_$level
DimPlot(Mac)
DimPlot(Mac, split.by = "orig.ident",ncol = 2,cols = c("grey","firebrick","black","royalblue"))

DotPlot(Mac, features = c("Itga4","P2ry12","Tmem119","H2-Aa","H2-DMb1","Klra2","Fxyd5","Maf","Mrc1","Cd163","Gas6","Ccrl2","Ccl3","Ccl4","Il1b","Nfkbia","Lgals3","Ccr2","Vim","Folr2","Ccl7","Ccr1","Lyve1","Ecrg4","Ttr","Enpp2")) + 
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "lightblue2",high= "red2",low="darkblue")

#Adding microglia and monocyte cluster for comparison as well
Myeloid.combined <- merge(Mac, y = MG, add.cell.ids = c("Mac", "MG"), project = "MyeloidCombined")
Myeloid.combined.2 <- merge(Myeloid.combined, y = Mo, add.cell.ids = c("Myeloid.combined", "Mo"), project = "MyeloidCombined")
levels(Myeloid.combined.2) = c("MG","Mo","MDM","MHCIIhi BAM","MHCIIlo BAM_1","MHCIIlo BAM_2")
DotPlot(Myeloid.combined.2, dot.scale = 9,features = c("Fxyd5","Lgals3","Ccr2","Vim","H2-Aa","Klra2","Pf4",
                                         "Ccrl2","Ccl3","Ccl4","Il1b","Mrc1","Cd163","Gas6","Folr2","P2rx7","Lilra5")) +
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "lightblue2",high= "red2",low="darkblue")


MHCIIhi <- WhichCells(Mac, idents = c("MHCIIhi BAM","MDM"))
MHCIIlow <- WhichCells(Mac, idents = c("MHCIIlo BAM_1","MHCIIlo BAM_2"))
DimPlot(Mac, label=F, cells.highlight= list(MHCIIlow, MHCIIhi), cols.highlight = c("cyan3","mediumorchid"))



###---------------------------------------------------------------###
#######---------Trajectory analysis running SCORPIUS--------------#########
###---------------------------------------------------------------###


Mac_2 <- subset(Mac, idents=c("MHCIIlo BAM_1","MHCIIlo BAM_2","MDM","MHCIIhi BAM"))
DimPlot(Mac_2)
library(SCORPIUS)
#Load normalized counts
DefaultAssay(Mac_2) <- "SCT"
Mac_2 <- NormalizeData(Mac_2)
expression <- t(Mac_2@assays$RNA@data)
expression_SCT <- t(Mac_2@assays$SCT@data)
group_name <- Mac_2@active.ident

#using norm RNA assay
space <- reduce_dimensionality(expression, dist = "spearman", ndim = 3)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space)

my_colour_palette <- c("MHCIIhigh BAM" = "grey","MDM" = "firebrick","MHCIIlo BAM_1" = "black", "MHCIIlo BAM_2" = "royalblue")
draw_trajectory_plot(
  space, 
  progression_group = group_name,progression_group_palette = my_colour_palette,
  path = traj$path,
  contour = TRUE, point_size = 1
)

#using SCT assay
space_sct <- reduce_dimensionality(expression_SCT, dist = "spearman", ndim = 3)
draw_trajectory_plot(space_sct, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space_sct)

#draw_trajectory_plot(space[, c(1, 3)]) + labs(y = "Component 3")
#draw_trajectory_plot(space[, c(1, 3)], progression_group = group_name) + labs(y = "Component 3")

#draw_trajectory_plot(space[, c(2, 3)]) + labs(x = "Component 2", y = "Component 3")
#draw_trajectory_plot(space[, c(2, 3)], progression_group = group_name) + labs(x = "Component 2", y = "Component 3")

expression <- as.matrix(expression_SCT)
gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:100,]
gene_sel <- gene_sel[!grepl("Rp", x=gene_sel$gene),]

expr_sel <- expression[,gene_sel$gene]
draw_trajectory_heatmap(expr_sel, traj$time, group_name, show_labels_row = T)
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules,show_labels_row = T,progression_group_palette = my_colour_palette)

),
                           minSamples=ncol(exprMat)*.01)

# (Adjust minimum values according to your dataset)


  ###---------------------------------------------------------------###
#######-----------------------DE Analysis---------------------#########
###---------------------------------------------------------------###

  Mac_MDM <- subset(Mac, idents=c("MDM"))
  Idents(Mac_MDM)=Mac_MDM$orig.ident
deg_Mac_MDM <-FindMarkers(object = Mac_MDM,
                      ident.1 = "ATBI_aCD49d",
                      ident.2 = "ATBI_Isotype",
                      test.use = "MAST",
                      logfc.threshold = -Inf,
                      min.pct = 0.5,
                      assay = "RNA"
) 

# Remove ribosomal, mitochondrial, and HLA genes
deg_Mac_MDM <- deg_Mac_MDM[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(deg_Mac_MDM)),]
deg_Mac_MDM <- deg_Mac_MDM[!grepl("AU", rownames(deg_Mac_MDM)),]
deg_Mac_MDM <- deg_Mac_MDM[!grepl("Gm", rownames(deg_Mac_MDM)),]

# Run Benjamini-Hochberg adjustment
deg_Mac_MDM$BH <- p.adjust(deg_Mac_MDM$p_val, method = "BH")
deg_Mac_MDM$genes <- rownames(deg_Mac_MDM)

# Create volcano plot
options(ggrepel.max.overlaps = 10)
volcano_plot(deg_Mac_MDM, title = "MDM with aCD49d Ab Treatment in Aged Mice PTBI",
             padj.thresh = 0.01, lfc.thresh = 0.25)+ scale_color_manual(values = c("red" = "blue"))+
#+ xlim(-0.5,2)
#+ ylim(0,5)
#Repeat for other subsets
