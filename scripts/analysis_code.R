# Derived from: https://swbioinf.github.io/spatialsnippets/d_cosmxIBD.html
library(Seurat)
library(SeuratObject)
library(tidyverse)
# requred: sf, harmony,  presto
# also:
# patchwork, DT, limma, edgeR
# others tba




## Paths
data_dir              <- file.path("data/")
raw_data_dir          <- file.path("raw_data/")

# Full size data
seurat_file_00_raw    <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_00_raw.RDS")

# Subsampled data
seurat_file_00_raw_subset          <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS")
seurat_file_01_preprocessed_subset <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS")
seurat_file_02_labelled_subset     <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_02_labelled_subsampled.RDS")

## Edited seurat load funcionts

# See comments here
#https://github.com/satijalab/seurat/discussions/9261
# base load nanostring
source("scripts/LoadNanostring_edited_function.R")



##  config

max_pc_negs        <- 1.5
max_avg_neg        <- 0.5

sample_codes <- c(HC="Healthy controls",UC="Ulcerative colitis",CD="Crohn's disease")


################################################################################

## Load one sample?
sample_path = file.path(raw_data_dir, "GSM7473682_HC_a")
so <- LoadNanostring(sample_path,
                     assay='RNA',
                     fov="GSM7473682_HC_a")

#NB: This default method drops most of the metatdata in the seurat object.
# e.g. what fov is each cell a member of? is missing.
# So we won't actually Use it.
# Trying to get this fixed in seurat.

#so@meta.data
#orig.ident nCount_RNA nFeature_RNA
#1_1   SeuratProject        368          189
#2_1   SeuratProject        810          286
#3_1   SeuratProject        119           74
so$tissue_sample   <- "HC_A"
so$group           <- "HC"
so$condition       <- "Healthy Controls"


################################################################################
# Load muliple samples
# Example code here, but do not run.


################################################################################
# Load a subsampled dataset for working with
# That takes too long, so loading a preloaded object

so.raw <- readRDS(seurat_file_00_raw_subset)

################################################################################
# Spatial plot
################################################################################

# Subset a seurat objecct by cells, just like a table.
# so[features, cells]
so.sample <- so.raw[, so.raw$tissue_sample=="HC_a"]

#so.sample
#An object of class Seurat
#999 features across 8795 samples within 2 assays
#Active assay: RNA (980 features, 0 variable features)
#1 layer present: counts.GSM7473682_HC_a
#1 other assay present: negprobes
#1 spatial field of view present: GSM7473682_HC_a

# seurat fov = slide
# bruker cosmx fov = region on slide

ImageDimPlot(so.sample,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             group.by = "fov_name",
             boundaries   = "segmentation",
             crop=TRUE,
             nmols = 10000)


# Also we have 'centroids' - whats the difference?
ImageDimPlot(so.sample,
             fov          = "GSM7473682_HC_a",
             axes = TRUE,
             border.color = "white", border.size = 0.1,
             group.by = "fov_name",
             boundaries   = "centroids",
             crop=TRUE,
             nmols = 10000)


#And and individual gene
ImageDimPlot(so.sample,
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
so.raw <- JoinLayers(so.raw)
GetAssayData(so.raw, assay = "RNA", layer = "counts")


################################################################################
# QC + Filtering
################################################################################


# Total counts per cell
ggplot(so.raw@meta.data, aes(x=nCount_RNA, col=orig.ident)) +
  geom_density() +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Transcript Counts per cell")

# Where might we put a threhsold?
min_count_per_cell <- 100



## Negative probe counts

#Negative counts per cell
#Negatvive probes are in a separate assay.
# This is a matter of preference, you could keep them with the rest.
so[['RNA']]
so[['negprobes']]

so.raw$pc_neg <-  ( so.raw$nCount_negprobes / (so.raw$nCount_RNA + so.raw$nCount_negprobes) ) * 100
so.raw[["negprobes"]] <- JoinLayers(so.raw[["negprobes"]]) # For caluclating these, need to have the negprobes merged
so.raw$avg_neg <-  colMeans(so.raw[["negprobes"]])   # only defined firsts sample.


# Discuss difference with Xenium - xenium lower counts, lower background.
max_pc_neg         <- 5

ggplot(so.raw@meta.data, aes(y=avg_neg, x=nCount_RNA)) +
  geom_point(pch=3, alpha=0.1) +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Negative probes vs counts")

ggplot(so.raw@meta.data, aes(y=pc_neg, x=nCount_RNA)) +
  geom_point(pch=3, alpha=0.1) +
  geom_hline(yintercept = max_pc_neg, lty=3) +
  geom_vline(xintercept = min_count_per_cell, lty=3) +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Negative probes vs counts")




### Apply a filter

# Should we apply a filter? High stringnec makes the analysis 'easier', but you get gaps spatially.

dim(so.raw)
so <- so.raw[ ,so.raw$nCount_RNA >= min_count_per_cell & so.raw$pc_neg <= max_pc_neg ]
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



