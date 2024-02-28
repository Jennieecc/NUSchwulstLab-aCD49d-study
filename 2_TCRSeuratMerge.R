# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
})

# Initialize paths
seurat_object <- "~"
output_dir <- "~"

# Source helper functions
source("~/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# General formatting fixes

# Load in seurat object
load(seurat_object)

# Add TCR Data
# Load in TCR information
contigs_merged <- read.csv("~contigs_merged.csv")
TCRs_paired <- read.csv("~paired_clonotypes_merged.csv")

# Remove duplicate rows per cell
tcr <- contigs_merged[!duplicated(contigs_merged$barcode), ]

# Only keep the barcode and clonotype columns
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode", "raw_clonotype_id","ID")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

# Slap the AA sequences onto our original table by clonotype_id
tcr <- merge(tcr, TCRs_paired[, c("clonotype_id", "trb_cdr3s", "tra_cdr3s", "frequency")])

# Reorder so barcodes are first column and set them as rownames.
tcr <- tcr[, c(2,1,3,4,5,6)]
rownames(tcr) <- tcr[,1]
tcr[,1] <- NULL
print("Head of TCR")
head(tcr)

# Add to the Seurat object's metadata.
samples_integrated_CD49d <- readRDS(file = "~samples_integrated_CD49d.rds")
samples_integrated_CD49d <- AddMetaData(object=samples_integrated_CD49d,
                 metadata=tcr)

# Add clonal info to metadata
samples_integrated_CD49d$clonal <- "NA"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency > 1] <- "C"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency == 1] <- "NC"
table(samples_integrated_CD49d@meta.data$clonal)

# Divide clonality by the frequency of clonotypes
samples_integrated_CD49d$clonal <- "NA"

samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency == 1] <- "Single"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency >1 &  samples_integrated_CD49d$frequency <5] <- "Small"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency >=5 &  samples_integrated_CD49d$frequency <40] <- "Medium"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency >=40 &  samples_integrated_CD49d$frequency <=101] <- "Large"
samples_integrated_CD49d$clonal[
  samples_integrated_CD49d$frequency > 101] <- "Hyperexpanded"

table(samples_integrated_CD49d@meta.data$clonal)



# Divide clonality by the frequency of clonotypes
query.projected_clonaltype$clonal <- "NA"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency == 1] <- "Single"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >1 &  query.projected_clonaltype$frequency <50] <- "Medium"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >=2 &  query.projected_clonaltype$frequency <35] <- "Medium"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency >=35 &  query.projected_clonaltype$frequency <=101] <- "Large"

query.projected_clonaltype$clonal[
  query.projected_clonaltype$frequency > 101] <- "Hyperexpanded"

table(query.projected_clonaltype@meta.data$clonal)

#------------------------------------------------------------------------------
# Filter TCRs

# Remove TCR annotations from Non T cells
t_cell_clusters <- c("T")

print(table(s[[c("cell_type_8_groups", "clonal")]]))

s@meta.data$clonal[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$clonotype_id[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$trb_cdr3s[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$tra_cdr3s[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$frequency[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$disease.clonal[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA
s@meta.data$Genderclonal[
  which(s@meta.data$cell_type_8_groups %!in% t_cell_clusters)] <- NA

print(table([[c("cell_type_8_groups", "clonal")]]))

# Calculate filtered frequency
filtered_freq <- data.frame(table(s[["clonotype_id"]]))
names(filtered_freq) <- c("clonotype_id", "frequency_filtered")
filtered_freq_merged <- plyr::join(s[["clonotype_id"]], filtered_freq)

# If in same order, merge
if (all.equal(filtered_freq_merged$clonotype_id, s[["clonotype_id"]]$clonotype_id)) {
  s@meta.data$frequency_filtered <- filtered_freq_merged$frequency_filtered
} else {
  stop("Unable to merge.")
}

# Reassign NC
s@meta.data$clonal[
  which(s@meta.data$frequency_filtered == 1)] <- "NC"
print("After NC reassignmnet")
print(table(s[[c("cell_type_8_groups", "clonal")]]))

#------------------------------------------------------------------------------
# Generate normalized frequency

# Calculate TCR counts per sample. s is the object after singleR annotation
tcr_count <- s[[c("orig.ident", "clonotype_id", "frequency_filtered")]] %>%
  dplyr::filter(!is.na(frequency_filtered)) %>%
  dplyr::distinct(clonotype_id, .keep_all = TRUE) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarize(tcr_count = sum(frequency_filtered)) %>%
  dplyr::distinct()
head(tcr_count)

# Calculate median TCR counts
median_tcr_count <- median(tcr_count$tcr_count)
median_tcr_count

# Initialize column
s@meta.data$normalized_frequency <- rep(NA, nrow(s@meta.data))

# Normalize frequency for each sample
for (sample in tcr_count$orig.ident) {
  sample_tcr_count <- tcr_count %>% dplyr::filter(orig.ident == sample) %>% dplyr::pull(tcr_count)
  s@meta.data$normalized_frequency[
    which(s@meta.data$orig.ident == sample) ] <- s@meta.data$frequency_filtered[
      which(s@meta.data$orig.ident == sample) ] * (median_tcr_count / sample_tcr_count)
}

# Check output
head(s[[c("frequency", "frequency_filtered", "normalized_frequency")]])

# Save modified seurat object
save(s, file = "~s_tcrclean_new.rds")

