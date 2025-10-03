# dev_setup

Instructions for setting up the workshop.

## Data

Download the .RDS files from the 
[Google Drive](https://drive.google.com/drive/u/1/folders/1eg1PmwGSWp_63p2IxZESmeq12U_GgjDJ):

| File                                                      | Size     |
| --------------------------------------------------------- | -------- |
| GSE234713_CosMx          | 538.8 MB |
| GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS          | 538.8 MB |
| GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS | 616.1 MB |
| GSM7473682_HC_a.tar.gz                                    | 344.7 MB |

```bash
git clone git@github.com:swbioinf/intro-spatial-transcriptomics-workshop.git
cd intro-spatial-transcriptomics-workshop
mkdir -p data
```

Move data files in correct directories.

```bash
mv GSE* data

tar -xvf GSM7473682_HC_a.tar.gz
mv GSM7473682_HC_a raw_data
rm GSM7473682_HC_a.tar.gz
```

The data structure should look like:

```bash
tree raw_data data
```

```console
raw_data
└── GSM7473682_HC_a
    ├── cell_boundaries_sf.parquet
    ├── GSM7473682_HC_a_exprMat_file.csv
    ├── GSM7473682_HC_a_fov_positions_file.csv
    ├── GSM7473682_HC_a_metadata_file.csv
    ├── GSM7473682_HC_a-polygons.csv
    ├── GSM7473682_HC_a_tx_file.csv
    └── GSM7473682_HC_a_tx_file.parquet
data
├── GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS
├── GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS
├── GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS
└── GSE234713_CosMx_IBD_seurat_02_filtered_subset.RDS

```

## Packages

In Rstudio, open project `intro-spatial-transcriptomics-workshop.Rproj`

```r
# Install renv 
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Update renv version
renv::record("renv@1.1.4") 
renv::init(bioconductor = T) 

# Manually download some dependencies
remotes::install_github("cellgeni/schard")
renv::install(c("Rfast2", "ape"))

#activate renv to load all other packages
renv::restore()
renv::activate()
```

May need to synchronise installed packages and lockfile packages - refer to
in-line renv suggestions to resolve. For example, if a package is used but
not recorded, run `renv::snapshot()` and select install packages -> update
snapshot.

## Rendering

```r
bookdown::render_book("index.Rmd", "bookdown::gitbook")
```

