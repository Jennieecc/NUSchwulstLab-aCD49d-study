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


Mac_ <-readRDS(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/Mac_.rds")

Mac <- subset(s, idents=c("MΦ"))
Mac <- RunPCA(object = Mac, verbose = FALSE)
Mac <- RunUMAP(Mac, dims = 1:30, verbose = FALSE)
Mac <- FindNeighbors(object = Mac, dims = 1:30)
Mac <- FindClusters(object = Mac, resolution = 0.8) 
DimPlot(Mac)
Mac.markers.RNA <- FindAllMarkers(Mac_,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top15_comb <- Mac.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)

Mac_ <- subset(Mac, idents=c("0","2","3","4"))
Mac_ <- RunPCA(object = Mac_, verbose = FALSE)
Mac_ <- RunUMAP(Mac_, dims = 1:30, verbose = FALSE)
Mac_ <- FindNeighbors(object = Mac_, dims = 1:30)
Mac_ <- FindClusters(object = Mac_, resolution = 1.2) 
DimPlot(Mac_)

background <- FindAllMarkers(samples_integrated_CD49d,assay="RNA",slot="data", min.pct = 0.25, logfc.threshold = 0)
write.csv(background, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/background.csv")

all.markers.RNA <- FindAllMarkers(Mac_,assay="RNA",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top15_comb <- all.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 15, avg_log2FC)
Mac_$level <- plyr::mapvalues(Idents(Mac_), from=c("0","1","2","3","4","5","6"),
                              to=c("1","MDM","2","Meningeal BAM","1","Meningeal BAM","Ependymal cells"))

Mac_$level <- plyr::mapvalues(Idents(Mac_), from=c("1","MDM","2","Meningeal BAM","Ependymal cells"),
                              to=c("MHCIIhi BAM","MDM","MHCIIlo BAM_1","MHCIIlo BAM_2","Ependymal cells"))
Idents(Mac_)=Mac_$level
Mac_2 <- subset(Mac_, idents=c("MHCIIhi BAM","MDM","MHCIIlo BAM_1","MHCIIlo BAM_2"))

ggsave("test.tiff", units="in", width=6, height=3.5, dpi=300, compression = 'lzw')

#using tiff() and dev.off
tiff("test.tiff", units="in", width=6, height=3.5, dpi=300,compression = 'lzw')

DimPlot(Mac_2, split.by = "orig.ident",ncol = 2,cols = c("grey","firebrick","black","royalblue"))

MDM <- subset(Mac_2, idents=c("MDM"))
BAM_2 <- subset(Mac_2, idents=c("MHCIIlo BAM_2"))
BAM_1 <- subset(Mac_2, idents=c("MHCIIlo BAM_1"))
BAM_hi <- subset(Mac_2, idents=c("MHCIIhi BAM"))

DotPlot(Myeloid.combined, features = c("Itga4","Vcam1"),split.by = "orig.ident",cols="RdBu",dot.scale=9,
        scale.min = 10,scale.max = 60,col.max = 1, col.min = -1)
DotPlot(Mac_2, features = c("Ccrl2","Ccl3","Ccl4","Il1b"),split.by = "orig.ident",cols="RdBu",
        scale.min = 10,scale.max = 60,col.max = 2, col.min = -0.5,dot.scale = 9)

DotPlot(Mac_, features = c("Itga4","P2ry12","Tmem119","H2-Aa","H2-DMb1","Klra2","Fxyd5","Maf","Mrc1","Cd163","Gas6","Ccrl2","Ccl3","Ccl4","Il1b","Nfkbia","Lgals3","Ccr2","Vim","Folr2","Ccl7","Ccr1","Lyve1","Ecrg4","Ttr","Enpp2")) + 
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "lightblue2",high= "red2",low="darkblue")

dev.off()

MG <- subset(s, idents=c("MG"))
MG <- RunPCA(object = MG, verbose = FALSE)
MG <- RunUMAP(MG, dims = 1:30, verbose = FALSE)
MG <- FindNeighbors(object = MG, dims = 1:30)
MG <- FindClusters(object = MG, resolution = 1) 

Idents(MG)=MG$cell_type_8_groups
Myeloid.combined <- merge(Mac_2, y = MG, add.cell.ids = c("Mac_", "MG"), project = "MyeloidCombined")
Myeloid.combined.2 <- merge(Myeloid.combined, y = MG, add.cell.ids = c("Myeloid.combined", "MG"), project = "MyeloidCombined")
levels(Myeloid.combined.2) = c("MG","Mo","MDM","MHCIIhi BAM","MHCIIlo BAM_1","MHCIIlo BAM_2")
DotPlot(Myeloid.combined.2, dot.scale = 9,features = c("Fxyd5","Lgals3","Ccr2","Vim","H2-Aa","Klra2","Pf4",
                                         "Ccrl2","Ccl3","Ccl4","Il1b","Mrc1","Cd163","Gas6","Folr2","P2rx7","Lilra5")) +
  RotatedAxis() +
  theme(legend.position="right")+
  scale_color_gradient2(mid = "lightblue2",high= "red2",low="darkblue")


pt <- table(Idents(Mac_2), Mac_2$orig.ident)
pt <- as.data.frame(pt)
pt_ATBI_aCD49d = pt[c(1:4),]
pt_ATBI_Isotype = pt[c(5:8),]
pt_YTBI_aCD49d = pt[c(9:12),]
pt_YTBI_Isotype = pt[c(13:16),]

pt_ATBI_aCD49d$Var1 <- as.character(pt_ATBI_aCD49d$Var1)
pt_ATBI_Isotype$Var1 <- as.character(pt_ATBI_Isotype$Var1)
pt_YTBI_aCD49d$Var1 <- as.character(pt_YTBI_aCD49d$Var1)
pt_YTBI_Isotype$Var1 <- as.character(pt_YTBI_Isotype$Var1)

library(RColorBrewer)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values = brewer.pal(10, "Paired"))


MHCIIhi <- WhichCells(Mac_, idents = c("MHCIIhi BAM","MDM"))
MHCIIlow <- WhichCells(Mac_, idents = c("MHCIIlo BAM_1","MHCIIlo BAM_2"))
DimPlot(Mac_, label=F, cells.highlight= list(MHCIIlow, MHCIIhi), cols.highlight = c("cyan3","mediumorchid"))
pdf(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/My Plot.pdf",   # The directory you want to save the file in
    width = 4.6, # The width of the plot in inches
    height = 3.3) # The height of the plot in inches



###---------------------------------------------------------------###
#######---------Trajectory analysis running SCORPIUS--------------#########
###---------------------------------------------------------------###


Mac_TA <- subset(Mac_, idents=c("1","2","MDM","Meningeal BAM"))
DimPlot(Mac_TA)
library(SCORPIUS)
#Load normalized counts
DefaultAssay(Mac_TA) <- "SCT"
Mac_TA <- NormalizeData(Mac_TA)
expression <- t(Mac_TA@assays$RNA@data)
expression_SCT <- t(Mac_TA@assays$SCT@data)
group_name <- Mac_TA@active.ident

#using norm RNA assay
space <- reduce_dimensionality(expression, dist = "spearman", ndim = 3)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space)

my_colour_palette <- c("1" = "grey","MDM" = "firebrick","2" = "black", "Meningeal BAM" = "royalblue")
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

my_colour_palette <- c("1" = "darkgrey","MDM" = "firebrick","2" = "black", "Meningeal BAM" = "royalblue")
df(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/My Plot1.pdf",   # The directory you want to save the file in
   width = 5, # The width of the plot in inches
   height = 3) # The height of the plot in inches
draw_trajectory_plot(
  space_sct, 
  progression_group = group_name,progression_group_palette = my_colour_palette,
  path = traj$path,
  contour = TRUE, point_size = 0.5, contour_alpha = 0.25, 
) 
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
pdf(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/My Plot11.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 10) # The height of the plot in inches
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules,show_labels_row = T,progression_group_palette = my_colour_palette)

###---------------------------------------------------------------###
#######-------------comparing pseudomodules-----------------#########
###---------------------------------------------------------------###
suppressMessages({
  library(SCPA)
  library(dyno)
  library(Seurat)
  library(tidyverse)
  library(magrittr)
  library(ComplexHeatmap)
  library(circlize)
})

###---------------------------------------------------------------###
#######------------------running SCENIC----------------------#########
###---------------------------------------------------------------###


suppressMessages({
  library("R.utils")
  library("utils")
  library("graphics")
  library("stats")
  library("data.table")
  library("mixtools")
  library("GSEABase")
  library("methods")
  library("Biobase")
  library("zoo")
  library("DT")
  library("NMF")
  library("plotly")
  #library("BiocStyle")
  library("rmarkdown")
  library("doMC")
  library("doRNG")
  library("doParallel")
  library("foreach")
  library("dynamicTreeCut")})
suppressMessages({
  library("GENIE3")
  library("AUCell")
  library("RcisTarget")
  library("RcisTarget.mm9.motifDatabases.20k")
  library("SCENIC")
})

#Firt time analysis download mm9 dataset
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")

#Loading Data from seurat object
exprMat <- GetAssayData(Mac_TA, assay = "RNA",slot = "counts")
exprMat <- as.matrix(exprMat)
#cell information
cellInfo <- data.frame(CellType=Idents(Mac_TA))
cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$CellType))
saveRDS(cellInfo, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/cellInfo.Rds")
colVars <- list(CellType=c("1"="grey", 
                           "MDM"="firebrick", 
                           "2"="black", 
                           "Meningeal BAM"="royalblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#Initialize SCENIC settings
library(SCENIC)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")
# mc9nr: Motif collection version 9: 24k motifs
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}


org <- "mgi" # mgi for mouse
dbDir <- "/Users/zhangyingchen" # RcisTarget databases location
myDatasetTitle <- "Mac" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
motifAnnotations_mgi <- motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
saveRDS(scenicOptions, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/scenicOptions.Rds")


genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)
corrMat <- cor(t(exprMat_filtered), method="spearman")
save(corrMat, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/corrMat.RData")


################################################################################
##########RUN SCENIC: Create co-expression modules (based on GENIE3 output)#####
################################################################################


linkList <- getLinkList(exprMat_filtered, threshold=0.001) # (slighly faster)
colnames(linkList) <- c("TF", "Target", "weight")

#Order by weight
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
save(linkList, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/GENIE3_linkList.RData")


########## Creating TF modules ##########
##### A. Keep only the TF-targets with weight > 0.001 #####
quantile(linkList$weight, probs=c(0.75, 0.90))


plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
     ylab="Weight", xlab="Links sorted decreasingly")
abline(h=0.001, col="blue") # Threshold

###How many percent of the weight matrix survives this treshold?
sum(linkList$weight>0.001)/nrow(linkList)*100
# 100

###Do the filtering
linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
# Number of links over the threshold: 
nrow(linkList_001) 
# 1572955 = all links

##### B. Create the geneSets for each TF #####
## For each TF, select:
##1.Targets with weight > 0.001
##2.Targets with weight > 0.005
##3.Top 50 targets (targets with highest weight)
##4.Targets for which the TF is within its top 5 regulators
##5.Targets for which the TF is within its top 10 regulators
##6.Targets for which the TF is within its top 50 regulators

tfModules <- list()
linkList_001$TF <- as.character(linkList_001$TF)
linkList_001$Target <- as.character(linkList_001$Target)
head(linkList_001)

#### Create TF-modules:
# 1: Weight > 0.001 (filtered in previous step) 
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))

# 3: Top 50 targets for each TF
# ("w001" should be ordered decreasingly by weight)
tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])

# 4-6: Top regulators per target 
# (linkList_001 should be ordered by weight!)
linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))
save(linkList_001_byTarget, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/linkList_001_byTarget.RData")

nTopTfs <- c(5, 10, 50)
nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

library(reshape2); library(data.table)
topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
  nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
  melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
})
topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
topTFsperTarget.asDf <-  data.frame(rbindlist(topTFsperTarget, idcol=TRUE))
head(topTFsperTarget.asDf)
colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")

# Merge the all the gene-sets:
tfModules.melted <- melt(tfModules)
colnames(tfModules.melted) <- c("Target", "TF", "method")
tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)

save(tfModules, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/tfModules.RData")

# Basic counts:
rbind(nGeneSets=nrow(tfModules), 
      nTFs=length(unique(tfModules$TF)), 
      nTargets=length(unique(tfModules$Target)))
# nGeneSets 3599860
# nTFs        8175
# nTargets     944

########## Split into positive- and negative-correlated targets ##########

### Keep only correlation between TFs and potential targets
tfs <- unique(tfModules$TF)
corrMat <- corrMat[tfs,]

### Split TF modules according to correlation
##if corr > 0.03 => get value of 1 (is activation)
##if corr < -0.03 => get value of -1 (is repression)
##otherwise => get value of 0 (not able to say if the link between TF and target is an activation or repression)
tfModules_byTF <- split(tfModules, factor(tfModules$TF))
tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
{
  tf <- unique(tfGeneSets$TF)
  targets <- tfGeneSets$Target
  cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
save(tfModules_withCorr, file="int/1.7_tfModules_withCorr.RData")

load("int/1.7_tfModules_withCorr.RData")
head(tfModules_withCorr)
#           Target            TF method corr
# 1           Cck 1810024B03Rik   w001    1
# 2          Iars 1810024B03Rik   w001    0
# 3           Axl 1810024B03Rik   w001    0
# 4         Ndor1 1810024B03Rik   w001    0
# 5 1700096K18Rik 1810024B03Rik   w001    0
# 6 1110051M20Rik 1810024B03Rik   w001    0
table(tfModules_withCorr$method)
# top10perTarget          top50 top50perTarget  top5perTarget           w001           w005 
# 73290          31450         366450          36645        3045479         106373
table(tfModules_withCorr$corr)
# -1       0       1 
# 21974 2488722 1148991

dim(tfModules_withCorr)
# 3659687       4


#DE Analysis
  Mac_MeningealBAM <- subset(Mac_, idents=c("Meningeal BAM")) #later rename as MHCIIhigh BAM_2
  Idents(Mac_MeningealBAM)=Mac_MeningealBAM$orig.ident
  
  deg_BAM <-FindMarkers(object = Mac_MeningealBAM,
                     ident.1 = "ATBI_aCD49d",
                     ident.2 = "ATBI_Isotype",
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.5,
                     assay = "RNA"
  ) 
  
  # Remove ribosomal, mitochondrial, and HLA genes
  deg_BAM <- deg_BAM[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(deg_BAM)),]
  deg_BAM <- deg_BAM[!grepl("AU", rownames(deg_BAM)),]
  deg_BAM <- deg_BAM[!grepl("Gm", rownames(deg_BAM)),]
  
  # Run Benjamini-Hochberg adjustment
  deg_BAM$BH <- p.adjust(deg_BAM$p_val, method = "BH")
  deg_BAM$genes <- rownames(deg_BAM)
  
  # Write out results
  write.csv(deg_BAM, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/deg_BAM.csv")
  
  ggsave("test.tiff", units="in", width=4, height=4, dpi=600, compression = 'lzw')
  
  #using tiff() and dev.off
  tiff('test.tiff', units="in", width=4, height=4, res=600, compression = 'lzw')
  
  # Create volcano plot
  options(ggrepel.max.overlaps = 10)
  #use 00_volcanoMerge_function 
  volcano_plot(deg_BAM, title = "BAMs with aCD49d Ab Treatment in Aged Mice PTBI",
               padj.thresh = 0.01, lfc.thresh = 0.25)+ xlim(-0.5,2)
    #+ ylim(0,5)

  
  dev.off()
  
  #-------------------------------------------------------Repeat----------------------------------------
  
  Mac_MDM <- subset(Mac_2, idents=c("MDM"))
  Idents(Mac_MDM)=Mac_MDM$orig.ident
  
  degs_wilcox_Young <-FindMarkers(object = Mac_MDM,
                     ident.1 = "YTBI_aCD49d",
                     ident.2 = "YTBI_Isotype",
                     test.use = "wilcox",
                     logfc.threshold = -Inf,
                     min.pct = 0.5,
                     assay = "RNA"
  ) 
  #4Xnumber of genes, use cells instead individual sample
  
  # pseudobulk the counts based on donor-condition-celltype
  #pseudo_Mac_MDM <- AggregateExpression(Mac_MDM, assays = "RNA", return.seurat = T)
  #Idents(pseudo_Mac_MDM) <- "orig.ident"
  
  #bulk.MDM.de <- FindMarkers(object = pseudo_Mac_MDM, 
                              #ident.1 = "YTBI-aCD49d", 
                              #ident.2 = "YTBI-Isotype",
                              #test.use = "DESeq2")
  
  
  # Remove ribosomal, mitochondrial, and HLA genes
degs_wilcox <- degs_wilcox[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(degs_wilcox)),]
degs_wilcox <- degs_wilcox[!grepl("AU", rownames(degs_wilcox)),]
degs_wilcox <- degs_wilcox[!grepl("Gm", rownames(degs_wilcox)),]
  
  # Run Benjamini-Hochberg adjustment
degs_wilcox$BH <- p.adjust(degs_wilcox$p_val, method = "BH")
degs_wilcox$genes <- rownames(degs_wilcox)

degs_wilcox_Young$BH <- p.adjust(degs_wilcox_Young$p_val, method = "BH")
degs_wilcox_Young$genes <- rownames(degs_wilcox_Young)
  # Write out results
  write.csv(degs, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/deg_MDM_Young.csv")
  
  # Create volcano plot
  options(ggrepel.max.overlaps = 10)
  #use 00_volcanoMerge_function 
  volcano_plot(degs_wilcox_Young, title = "MDMs with aCD49d Ab Treatment in Aged Mice PTBI",
               padj.thresh = 0.01, lfc.thresh = 0.25)
    #ylim(0,5)
  
  volcano_plot(degs_wilcox_Young, title = "MDMs with aCD49d Ab Treatment in Aged Mice PTBI",
               padj.thresh = 0.01, lfc.thresh = 0.25)
  
  pdf(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/My Plot.pdf",   # The directory you want to save the file in
      width = 7, # The width of the plot in inches
      height = 7) # The height of the plot in inches
 
Idents(Mac_)=Mac_$orig.ident
Mac_Aged <- subset(Mac_, idents=c("ATBI_aCD49d","ATBI_Isotype"))
Idents(Mac_Aged)=Mac_Aged$level

Mac_Aged <- subset(Mac_Aged, idents=c("MDM","Meningeal BAM"))
Idents(Mac_Aged)=Mac_Aged$level
DoMultiBarHeatmap(Mac_Aged, features = c("Cd86","Ctsb","Pla2g7","Slc11a1","Malat1","Clec4n","Maf","Cd163","Clec12a","Klra2","Fxyd5","Ccr2","Lgals3",
                                     "Actg1","Fau","Tpt1"), additional.group.by = "orig.ident",size = 3,additional.group.sort.by = "orig.ident")+scale_fill_viridis(option = "J")
  

pt <- table(Idents(Mac_), Mac_$orig.ident)
pt <- as.data.frame(pt)
pt1 = pt[c(1:4,6:9,11:14,16:19),]

pt_ATBI_aCD49d = pt[c(1:4),]
pt_ATBI_Isotype = pt[c(6:9),]
pt_YTBI_aCD49d = pt[c(11:14),]
pt_YTBI_Isotype = pt[c(16:19),]

pt_ATBI_aCD49d$Var1 <- as.character(pt_ATBI_aCD49d$Var1)
pt_ATBI_Isotype$Var1 <- as.character(pt_ATBI_Isotype$Var1)
pt_YTBI_aCD49d$Var1 <- as.character(pt_YTBI_aCD49d$Var1)
pt_YTBI_Isotype$Var1 <- as.character(pt_YTBI_Isotype$Var1)


library(RColorBrewer)
ggplot(pt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values = brewer.pal(10, "Paired"))  

#pathway analysis result plot
ggsave("test.tiff", units="in", width=13, height=3.5, dpi=300, compression = 'lzw')

#using tiff() and dev.off
tiff("test.tiff", units="in", width=13, height=3.5, dpi=300, compression = 'lzw')

GO_Selected_Mac <- read.table(file = "/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/GO_selected.csv", 
                              sep = ",", header=TRUE)
ggplot(GO_Selected_Mac,aes(reorder(Description, Enrichment), Enrichment, fill=FDR.q.value, split='Description')) + geom_col(width = 0.1) + 
  theme(axis.text.x=element_text(angle=-40, hjust=0, size = 15), axis.text.y = element_text(size = 16))+
  RotatedAxis()+coord_flip(xlim=c(length(unique(GO_Selected_Mac$Description))-7,length(unique(GO_Selected_Mac$Description))))+
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 20))+geom_bar(stat = "identity",color="black") +scale_fill_distiller(palette = "Blues")+
  xlab("Top GO Terms Macrophages")


#-------------------------------------------------------remove ependemal cells and rename----------------------------------------

#DE Analysis
BAM_2 <- subset(Mac_2, idents=c("MHCIIlo BAM_2"))
Idents(BAM_2)=BAM_2$orig.ident

deg_lowBAM_2 <-FindMarkers(object = BAM_2,
                      ident.1 = "ATBI_aCD49d",
                      ident.2 = "ATBI_Isotype",
                      test.use = "MAST",
                      logfc.threshold = -Inf,
                      min.pct = 0.5,
                      assay = "RNA"
) 

# Remove ribosomal, mitochondrial, and HLA genes
deg_lowBAM_2 <- deg_lowBAM_2[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(deg_lowBAM_2)),]
deg_lowBAM_2 <- deg_lowBAM_2[!grepl("AU", rownames(deg_lowBAM_2)),]
deg_lowBAM_2 <- deg_lowBAM_2[!grepl("Gm", rownames(deg_lowBAM_2)),]

# Run Benjamini-Hochberg adjustment
deg_lowBAM_2$BH <- p.adjust(deg_lowBAM_2$p_val, method = "BH")
deg_lowBAM_2$genes <- rownames(deg_lowBAM_2)

# only Atf3, Nuak2, and Zfp36 were significantly upregulated no downregulated genes

# Create volcano plot
options(ggrepel.max.overlaps = 10)
#use 00_volcanoMerge_function 
volcano_plot(deg_lowBAM_2, title = "MHCII low BAMs with aCD49d Ab Treatment in Aged Mice PTBI",
             padj.thresh = 0.01, lfc.thresh = 0.25)+ scale_color_manual(values = c("red" = "blue"))+
 xlim(-0.5,2)
#+ ylim(0,5)
#-------------------------------------------------------repeat----------------------------------------

Mac_MHCIIhigh_BAM <- subset(Mac_2, idents=c("MHCIIhi BAM"))
Idents(Mac_MHCIIhigh_BAM)=Mac_MHCIIhigh_BAM$orig.ident

deg_highBAM <-FindMarkers(object = Mac_MHCIIhigh_BAM,
                         ident.1 = "ATBI_aCD49d",
                         ident.2 = "ATBI_Isotype",
                         test.use = "MAST",
                         logfc.threshold = -Inf,
                         min.pct = 0.5,
                         assay = "RNA"
) 

# Remove ribosomal, mitochondrial, and HLA genes
deg_highBAM <- deg_highBAM[-grep(pattern = "^Rps|^Rpl|^Mt-|^mt-|^Hla-|^Gm-|^rp-|^H2-", x = rownames(deg_highBAM)),]
deg_highBAM <- deg_highBAM[!grepl("AU", rownames(deg_highBAM)),]
deg_highBAM <- deg_highBAM[!grepl("Gm","Junb","Egr1","Fos", rownames(deg_highBAM)),]

# Run Benjamini-Hochberg adjustment
deg_highBAM$BH <- p.adjust(deg_highBAM$p_val, method = "BH")
deg_highBAM$genes <- rownames(deg_highBAM)

write.csv(deg_highBAM, file="/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/deg_highBAM_Aged.csv")

# Create volcano plot
options(ggrepel.max.overlaps = 10)

volcano_plot(deg_highBAM, title = "MHCII high BAMs with aCD49d Ab Treatment in Aged Mice PTBI",
             padj.thresh = 0.01, lfc.thresh = 0.25)+ xlim(-1,3)+ scale_color_manual(values = c("red" = "forestgreen"))
+ xlim(-0.5,3)

deg_merge_Mac <- read.csv("/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/_BAM_Merge_ATBIaCD49d_vs_ATBI_Isotype_degs.csv")
deg_merge_Mac <- read.csv("/Users/zhangyingchen/Desktop/AvY_antiCD49d_study/scRNA/_BAM_2_cell_Merge_ATBIaCD49d_vs_ATBI_Isotype_degs.csv")
#use 00_volcanoMerge_function 
volcano_plot(deg_merge_Mac, title = "CD8+ T Cells with aCD49d Ab 
  Treatment in Aged Mice PTBI",
             padj.thresh = 0.01, lfc.thresh = 1)
