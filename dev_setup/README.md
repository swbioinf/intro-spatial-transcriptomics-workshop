# dev_setup

Instructions for setting up the workshop.

## Data

| File                                                      | Size     |
| --------------------------------------------------------- | -------- |
| GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS          | 538.8 MB |
| GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS | 616.1 MB |
| GSM7473682_HC_a.tar.gz                                    | 344.7 MB |

## Packages

```r
# Install renv 
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Update renv version
renv::record("renv@1.1.4") 

# activate renv
renv::restore()
renv::activate()
```

May need to synchronise installed packages and lockfile packages - refer to
in-line renv suggestions to resolve.


