# Derived from: https://swbioinf.github.io/spatialsnippets/d_cosmxIBD.html
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)

# Full size data
seurat_file_00_raw    <- here("data", "GSE234713_CosMx_IBD_seurat_00_raw.RDS")
seurat_file_01_loaded <- here("data", "GSE234713_CosMx_IBD_seurat_01_loaded.RDS")

# Subsampled data
seurat_file_00_raw_subset <-  here("data", "GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS")





################################################################################

## Subsample raw data
# Take up to 4 FOVS from each sample in the HC and CD groups (6 samples)
if (! file.exists(seurat_file_00_raw_subset)) {

  so.raw.full <- readRDS(seurat_file_00_raw)

  # QUick and dirty - Seurat does not support "_" in fov names.
  # And by default will use the project sample names.
  # It will work, and do many things, but molecules won't plot.
  # Go through and reaname all of these
  for (fov_name in names(so.raw.full@images)) {
    # rename.
    fov_obj <- so.raw.full[[fov_name]]
    so.raw.full[[fov_name]] <- NULL
    so.raw.full[[gsub("_","-", fov_name)]] <- fov_obj

  }
  names(so.raw.full@images)


  so.raw.subset <- so.raw.full [, so.raw.full$group %in% c("CD", "HC") & so.raw.full@meta.data$fov <= 4]
  saveRDS(  so.raw.subset, seurat_file_00_raw_subset)

}






