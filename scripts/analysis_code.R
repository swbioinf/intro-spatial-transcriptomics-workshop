# Derived from: https://swbioinf.github.io/spatialsnippets/d_cosmxIBD.html
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
library(sf)
library(harmony)
library(remotes) # to install presto
# remotes::install_github("immunogenomics/presto")

# required:
# patchwork, DT, limma, edgeR
# others tba

# Full size data
seurat_file_00_raw    <- here("data", "GSE234713_CosMx_IBD_seurat_00_raw.RDS")

# Subsampled data
seurat_file_00_raw_subset          <- here("data", "GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS")
seurat_file_01_preprocessed_subset <- here("data", "GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS")
#seurat_file_02_labelled_subset     <- here("data", "GSE234713_CosMx_IBD_seurat_02_labelled_subsampled.RDS")

## Edited seurat load funcionts

# See comments here
#https://github.com/satijalab/seurat/discussions/9261
# base load nanostring
source("scripts/LoadNanostring_edited_function.R")



##  config
# TODO: explain choices
max_pc_negs        <- 1.5
max_avg_neg        <- 0.5

sample_codes <- c(HC="Healthy controls",UC="Ulcerative colitis",CD="Crohn's disease")


################################################################################

## Load one sample?
sample_path = here("raw_data", "GSM7473682_HC_a")
so <- LoadNanostring(sample_path,
                     assay='RNA',
                     fov="GSM7473682_HC_a")

# Explain input data, what it represents, how it will be used for analysis.
# Define terms like centroids, segments, fov

#NB: This default method drops most of the metadata in the seurat object.
# e.g. what fov is each cell a member of? is missing.
# So we won't actually Use it.
# Trying to get this fixed in seurat.

# Inspect Seurat object
# Explain that Seurat objects are a data structure to store the count data,
# and additional metadata. As the data is very high dimensional, it is important
# that connected information about cells, samples, and data stay connected
# during data analysis
so

#An object of class Seurat
#999 features across 42045 samples within 1 assay
#Active assay: RNA (999 features, 0 variable features)
#1 layer present: counts
#1 spatial field of view present: GSM7473682_HC_a

# Explain: what does each thing define? e.g. features, layers, fovs

# Inspect metadata
head(so@meta.data)
#orig.ident nCount_RNA nFeature_RNA
#1_1   SeuratProject        368          189
#2_1   SeuratProject        810          286
#3_1   SeuratProject        119           74

# Add metadata. Needed to compare
so$condition       <- "Healthy Controls"
so$tissue_sample   <- "HC_A" # denotes healthy control, from sample A

################################################################################
# Load muliple samples
# Example code here, but do not run.


################################################################################
# Load a subsampled dataset for working with multiple samples
# That takes too long, so loading a preloaded object
# Add figure that represents the subset

# Explain: Subset data, how many samples and conditions. Represent graphically
# e.g. three healthy controls (HC) vs. three CD
so_raw <- readRDS(seurat_file_00_raw_subset)
so_raw

# An object of class Seurat
# 999 features across 68624 samples within 2 assays
# Active assay: RNA (980 features, 0 variable features)
# 6 layers present: counts.GSM7473682_HC_a, counts.GSM7473683_HC_b, counts.GSM7473684_HC_c, counts.GSM7473688_CD_a, counts.GSM7473689_CD_b, counts.GSM7473690_CD_c
# 1 other assay present: negprobes
# 6 spatial fields of view present: GSM7473682_HC_a GSM7473683_HC_b GSM7473684_HC_c GSM7473688_CD_a GSM7473689_CD_b GSM7473690_CD_c
# Explain: difference between this SO vs. previous one (`so`)
#     - layers represent the different samples and their condition
#     - RNA and negprobes assays and what data they represent

# Display metadata columns
names(so_raw@meta.data)

# There is a lot of information. Use something like skimr for EDA and to
# understand each field
skimr::skim(so_raw@meta.data)

# Explain: number of rows (68624) represent data points (?)
# Identify important metadata columns and discuss what they mean

################################################################################
# Spatial plot
################################################################################

# View the tissue_sample metadata to prepare for the next step
table(so_raw$tissue_sample)

# Subset a seurat object by cells, just like a table.
# so[features, cells]
so_sample <- so_raw[, so_raw$tissue_sample=="HC_a"]

# Some warnings with not validating FOV, centroid, and seurat objects

so_sample
#An object of class Seurat
#999 features across 8795 samples within 2 assays
#Active assay: RNA (980 features, 0 variable features)
#1 layer present: counts.GSM7473682_HC_a
#1 other assay present: negprobes
#1 spatial field of view present: GSM7473682_HC_a

# seurat fov = slide
# bruker cosmx fov = region on slide
# Note: layer and FOVs are now only for a single sample, HC_a

# Visualise clusters in a spatial context
# Explain: key arguments
ImageDimPlot(so_sample,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             group.by = "fov_name",
             boundaries   = "segmentation",
             crop=TRUE,
             nmols = 10000)


# Also we have 'centroids' - whats the difference?
ImageDimPlot(so_sample,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             group.by = "fov_name",
             boundaries   = "centroids",
             crop=TRUE,
             nmols = 10000)


#And and individual gene
ImageDimPlot(so_sample,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             #group.by = "fov_name",
             boundaries   = "segmentation",
             molecules = "S100A6",
             crop=TRUE,
             nmols = 10000)

###################
# Where are my counts?
# GetAssayData(so.raw, assay = "RNA", layer = "counts")
# Oh no an error:
#Error in `GetAssayData()` at SeuratObject/R/seurat.R:1901:3:
#  ! GetAssayData doesn't work for multiple layers in v5 assay.



# Seurat *5* vs Seurat 4. Multiple assays fovs.
# SOme operations only work on a single
# Join layers / split layers.
# tHis is useful for later steps, but annoying now.
so_raw <- JoinLayers(so_raw)
GetAssayData(so_raw, assay = "RNA", layer = "counts")


################################################################################
# QC + Filtering
################################################################################

# Total counts per cell
ggplot(so_raw@meta.data, aes(x=nCount_RNA, col=orig.ident)) +
  geom_density() +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Transcript Counts per cell")

# Where might we put a threhsold?
min_count_per_cell <- 100



## Negative probe counts

#Negative counts per cell
#Negative probes are in a separate assay.
# This is a matter of preference, you could keep them with the rest.
so_raw[['RNA']]
so_raw[['negprobes']]

so_raw$pc_neg <-  ( so_raw$nCount_negprobes / (so_raw$nCount_RNA + so_raw$nCount_negprobes) ) * 100
so_raw[["negprobes"]] <- JoinLayers(so_raw[["negprobes"]]) # For caluclating these, need to have the negprobes merged
so_raw$avg_neg <-  colMeans(so_raw[["negprobes"]])   # only defined firsts sample.


# Discuss difference with Xenium - xenium lower counts, lower background.
max_pc_neg         <- 5

ggplot(so_raw@meta.data, aes(y=avg_neg, x=nCount_RNA)) +
  geom_point(pch=3, alpha=0.1) +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Negative probes vs counts")

ggplot(so_raw@meta.data, aes(y=pc_neg, x=nCount_RNA)) +
  geom_hline(yintercept = max_pc_neg, lty=3, colour = "red") +
  geom_vline(xintercept = min_count_per_cell, lty=3, colour = "red") +
  geom_point(pch=3, alpha=0.1) +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Negative probes vs counts")




### Apply a filter

# Should we apply a filter? High stringency makes the analysis 'easier', but you get gaps spatially.

dim(so_raw)
so <- so_raw[ ,so_raw$nCount_RNA >= min_count_per_cell & so_raw$pc_neg <= max_pc_neg ]
dim(so)
table(so@meta.data$orig.ident)

# Basic preprocessing
# Split layers out again
# https://satijalab.org/seurat/articles/seurat5_integration
so <- split(so, f = so$orig.ident)

# Run through preprocessing
so <- NormalizeData(so)

## Do per sample to mimic paper approach somewhat.
so <- FindVariableFeatures(so, nfeatures = 200)
VariableFeaturePlot(so)

# Explain: Each point represents a gene, red = HVGs selected

so <- ScaleData(so) # Just 200 variable features
so <- RunPCA(so, features = VariableFeatures(so))


DimPlot(so, reduction = "pca", group.by='tissue_sample')
DimPlot(so, reduction = "pca", group.by='tissue_sample', dims=c(3,4))

ElbowPlot(so)
num_dims <- 15


### No integrations => batch effect
demo_no_correction = FALSE
if (demo_no_correction)  {

  so <- RunUMAP(so, dims=1:num_dims, )

  DimPlot(so, group.by='tissue_sample') # If umap exists, umap.

  so <- FindNeighbors(so, dims = 1:num_dims)
  so <- FindClusters(so)   # takes a few min.
  # don't talk about optimisation of clusters, point people at single cell resources.

  DimPlot(so, group.by='seurat_clusters')
}

# Demo - supply people with output so
actually_run_harmony = FALSE
if (actually_run_harmony) {
  # https://satijalab.org/seurat/articles/seurat5_integration
  # Slow?
  so <- IntegrateLayers(so, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
  so <- RunUMAP(so, dims=1:num_dims, reduction = 'harmony')

  DimPlot(so, group.by='tissue_sample')

  so <- FindNeighbors(so, reduction = "harmony", dims = 1:num_dims)
  so <- FindClusters(so, resolution = 0.2)
  # don't talk about optimisation of clusters, point people at single cell resources.
}



# Load previous (or save)
if (! file.exists(seurat_file_01_preprocessed_subset)) {
  so <- readRDS(seurat_file_01_preprocessed_subset)
} else {
  saveRDS(so, seurat_file_01_preprocessed_subset)
}


DimPlot(so, group.by='seurat_clusters')
DimPlot(so, group.by='orig.ident')

# Cell counts per cluster per sample
table(so$orig.ident, so$seurat_clusters)
# < provide So object at this point>






################################################################################
## More plotting

DimPlot(so, group.by='orig.ident', split.by = 'condition')

so.sample <- so[, so$tissue_sample=="HC_a"]

table(so$fov_name)
so.fov    <- so[, so$fov_name=="HC_a_004"]

ImageDimPlot(so.fov,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             group.by = "seurat_clusters",
             boundaries   = "segmentation",
             cols="glasbey", #more distinct than default ggplot
             crop=TRUE,
             nmols = 10000)




## Marker genes
so <- JoinLayers(so)

Idents(so) <- so$seurat_clusters
marker_table <- FindAllMarkers(so, only.pos = TRUE, max.cells.per.ident = 100)

marker_table.top <- marker_table %>%
  filter(p_val_adj < 0.01)%>%
  select(gene, cluster, everything()) %>% # change column order
  group_by(cluster) %>%
  slice_min(order_by=p_val_adj, n=5)

marker_table.top


# List of the top gene for each cluster
marker_table.top_genes <- marker_table %>%
  filter(p_val_adj < 0.01)%>%
  select(gene, cluster, everything()) %>% # change column order
  group_by(cluster) %>%
  slice_min(order_by=p_val_adj, n=1) %>%
  pull(gene) %>%
  unique()

FeaturePlot(so, marker_table.top_genes, raster = TRUE)

VlnPlot(so, marker_table.top_genes, pt.size=0)
VlnPlot(so, 'nCount_RNA', pt.size=0) + scale_y_log10()



## Apply some cluster names

so$cluster_code <- factor( paste0("c", so$seurat_clusters),   levels=paste0('c', levels(so$seurat_clusters)))
Idents(so) <- so$cluster_code

# Will pick some sensible names later, onces settled.
# PIGR is a marker of intestinal epithelial, and very obviously epithelial.
cluster_content <- list(
  c0 = "Stromal",
  c1 = "IntestinalEpithelia",
  c2 = "CollagenRich",
  c3 = "Immuneish",
  c4 = "test",
  c5 = "fixme"
)

# c5 => c5: Mono
so$celltype <- factor (
  paste0(names(cluster_content[as.character(so$cluster_code)]), ": ", cluster_content[as.character(so$cluster_code)]) ,
  levels = paste0( names(cluster_content), ": ", cluster_content)
)



## Celltype proportions
celltype_summary_table<- so@meta.data %>% as_tibble() %>%
  group_by(condition, tissue_sample, celltype ) %>%
  summarise(n_cells = n())



ggplot(celltype_summary_table, aes(x=tissue_sample, y=n_cells, fill=celltype)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle( "Celltype composition") +
  facet_wrap(~condition, ncol = 1, scales = 'free_y')



# Next Steps:
# Differntial expression



