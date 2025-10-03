library(here)
library(SingleR)
#cellgeni/schard
library(schard)
library(SingleCellExperiment)
library(schard)
library(BiocParallel)
library(Seurat)

ts_largeintestine_h5ad <- file.path("data/tabula_sapiens/tabula_sapiens_large_intestine_82e3b450-6704-43de-8036-af1838daa7df.h5ad")
sce.ts_intestine = schard::h5ad2sce(ts_largeintestine_h5ad)

# Read existing data
so <- readRDS("data/GSE234713_CosMx_IBD_seurat_01_preprocessed_subsampled.RDS")
so <- JoinLayers(so)



sce.ts_intestine
table(sce.ts_intestine$cell_type)
table(sce.ts_intestine$broad_cell_class)
table(sce.ts_intestine$donor_id, sce.ts_intestine$broad_cell_class )


# Not needed, but first filter down to matched genes in our panel
sce.ts_intestine.genename <- sce.ts_intestine[rowData(sce.ts_intestine)$feature_name %in% rownames(so),]

# are there any duplicates (we'd need to handle them)
sum(table(rowData(sce.ts_intestine.genename)$feature_name) != 1 )

# just rename the genes to the gene names
rownames(sce.ts_intestine.genename) <-  rowData(sce.ts_intestine.genename)$feature_name

# Subsample
#cells2keep <- sample(size=1000, ncol(sce.ts_intestine.genename))
#sce.ts_intestine.genename <- sce.ts_intestine.genename[,cells2keep]




# Pull out the normalised matrix.
# Quirk of this coming from the python world, there's no counts assay, and the normalised is in 'X'
ref_norm_matrix <- assay(sce.ts_intestine.genename, 'X')


# Show a heatmap of htose at the sample summary level?
norm_matrix <- GetAssayData(so, assay = 'RNA', layer = 'data')

predictions_file <- file.path("data","predictions_with_ts_intestine_broad_cell_class.RDS")
predictions <- SingleR::SingleR(test = norm_matrix,
                                ref   = ref_norm_matrix,
                                labels = sce.ts_intestine.genename$broad_cell_class,
                                aggr.ref = TRUE, # builds a pseudubulk reference , speedier processing
                                BPPARAM = MulticoreParam(workers=8)
                                )
saveRDS(predictions, predictions_file)

print(" and detailed.")

predictions_file2 <- file.path("data","predictions_with_ts_intestine_detailed_celltype.RDS")
predictions2 <- SingleR::SingleR(test = norm_matrix,
                                ref   = ref_norm_matrix,
                                labels = sce.ts_intestine.genename$cell_type,
                                aggr.ref = TRUE, # builds a pseudubulk reference , speedier processing
                                BPPARAM = MulticoreParam(workers=8)
)
saveRDS(predictions2, predictions_file2)


#so$singleR_labels <- predictions$labels
#DimPlot(so, group.by='singleR_labels')









