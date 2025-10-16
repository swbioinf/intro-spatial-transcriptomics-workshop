# dev_setup

Instructions for setting up the workshop.

## Data

Files required for this workshop are stored in the developer
[google drive](https://drive.google.com/drive/u/1/folders/1eg1PmwGSWp_63p2IxZESmeq12U_GgjDJ).
These are used throughout the different .Rmd file (lessons) and are required
for the webpage to render.

## Directory setup

Retrieve a copy of the repository and set up `data` and `raw_data` folders.

```bash
git clone git@github.com:swbioinf/intro-spatial-transcriptomics-workshop.git && \
cd intro-spatial-transcriptomics-workshop && \
mkdir -p data raw_data
```

Download all data files from the
[google drive](https://drive.google.com/drive/u/1/folders/1eg1PmwGSWp_63p2IxZESmeq12U_GgjDJ).

Move data files in correct directories.

```bash
mv GSE* data/ && \
mv predictions* data/ && \
mv tabula_sapiens data/ && \
tar -xvf GSM7473682_HC_a.tar.gz && \
mv GSM7473682_HC_a/ raw_data/ && \
rm GSM7473682_HC_a.tar.gz
```

The data structure should look like:

```bash
tree raw_data data
```

```console
# Note: predictions files not required in this current version
data
├── GSE234713_CosMx_annotation.csv.gz
├── GSE234713_CosMx_IBD_seurat_00_raw_subsampled.RDS
├── GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS
├── GSE234713_CosMx_IBD_seurat_02_filtered_subset.RDS
├── GSE234713_CosMx_IBD_seurat_02_rna70_neg4.RDS
└── tabula_sapiens
    └── tabula_sapiens_large_intestine_82e3b450-6704-43de-8036-af1838daa7df.h5ad
raw_data
└── GSM7473682_HC_a
    ├── cell_boundaries_sf.parquet
    ├── GSM7473682_HC_a_exprMat_file.csv
    ├── GSM7473682_HC_a_fov_positions_file.csv
    ├── GSM7473682_HC_a_metadata_file.csv
    ├── GSM7473682_HC_a-polygons.csv
    ├── GSM7473682_HC_a_tx_file.csv
    └── GSM7473682_HC_a_tx_file.parquet
```

## Packages

TODO: Instructions will be updated for VM setup as the dependencies
and `renv.lockfile` are changing quickly.

In Rstudio, open project `intro-spatial-transcriptomics-workshop.Rproj`

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
in-line renv suggestions to resolve. For example, if a package is used but
not recorded, run `renv::snapshot()` and select install packages -> update
snapshot.

## Rendering 

Following best practices in other static site generator frameworks
(quarto, mkdocs), rendered files for the webpage will be kept in a separate
`gh-pages` branch. This focuses `main` and development  branches for source code
for ease of review, tracking, and minimising file conflicts.

Render the site locally. The following command stitches the separate .Rmd
files together, and generates the output .html files. These are
self-contained in the `docs/` directory.

```r
# Alternatively, in Rstudio: Build -> Build Book -> bookdown:bs4book
rmarkdown::render_site(output_format = 'bookdown::bs4_book', encoding = 'UTF-8')
```

The `docs/` directory, along with some intermediate files will be ignored by
git (`.gitignore`) to ensure they are not committed and pushed on the main
branch. When running `git status`, these should not appear as modified.

## Local testing

To review your changes in the rendered web page during development, open
`docs/index.html` in your browser.

## Publishing the page

### Configuring github pages

The repository should be configured so the `gh-pages` branch is used to publish the github pages.

Settings -> Pages -> Build and deployment:

- Source: Deploy from a branch
- Branch: `gh-pages`
- Folder: `/ (root)`
- Save

The next time a commit is pushed to the `gh-pages` branch, the updates
will be published.

### Setting up the worktree (first time only)

Git normally tracks one branch per folder. A worktree lets you have another
checkout (here, `gh-pages`) in a separate directory while keeping the main
branch clean.

```bash
git fetch origin gh-pages # get latest gh-pages branch
git worktree add ../rendered_docs gh-pages
```

Then verify.

```bash
git worktree list
```

The expected output should contain:

- The main repo
- A folder `rendered_docs/` containing the live `gh-pages` branch

e.g. 

```console
/home/fred/GitHub/intro-spatial-transcriptomics-workshop  <commit-hash> [main]
/home/fred/GitHub/rendered_docs                           <commit-hash> [gh-pages]
```

### Updating rendered files

To update the github page after building the book via Rstudio or `render_site()`
locally, replace the contents of the `../rendered_docs/` folder that tracks
the `gh-pages` branch.

```bash
rm -rf ../rendered_docs/*
cp -rv docs/* ../rendered_docs
```

Then, commit and push!

```bash
cd ../rendered_docs
git add --all
git commit -m "Update site"
git push gh-pages
```

## Converting content to working .Rmd files

Remove all text from an .Rmd file.

Example, for a single file.

```bash
scripts/strip_rmd.sh 45-normalisation.Rmd
```
```console
45-normalisation.Rmd → notebooks/45-normalisation.Rmd
```

For all .Rmd in the repo.

```bash
for i in *.Rmd; do scripts/strip_rmd.sh $i; done
```
```console
01-intro.Rmd → notebooks/01-intro.Rmd
05-background.Rmd → notebooks/05-background.Rmd
20-study.Rmd → notebooks/20-study.Rmd
22-example_output.Rmd → notebooks/22-example_output.Rmd
30-data.Rmd → notebooks/30-data.Rmd
40-QC.Rmd → notebooks/40-QC.Rmd
45-normalisation.Rmd → notebooks/45-normalisation.Rmd
50-ReduceDims.Rmd → notebooks/50-ReduceDims.Rmd
55-BatchCorrection.Rmd → notebooks/55-BatchCorrection.Rmd
60-clustering.Rmd → notebooks/60-clustering.Rmd
65-ClusterLabelling.Rmd → notebooks/65-ClusterLabelling.Rmd
72-spatiallyrestrictedgenes.Rmd → notebooks/72-spatiallyrestrictedgenes.Rmd
80-nextsteps.Rmd → notebooks/80-nextsteps.Rmd
85-resources.Rmd → notebooks/85-resources.Rmd
90-references.Rmd → notebooks/90-references.Rmd
index.Rmd → notebooks/index.Rmd
```
