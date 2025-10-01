
library(here)
library(SingleR)
#cellgeni/schard
library(schard)
library(SingleCellExperiment)


ts_largeintestine_h5ad <- file.path("data/tabula_sapiens/tabula_sapiens_large_intestine_82e3b450-6704-43de-8036-af1838daa7df.h5ad")


#
sce.ts_intestine = schard::h5ad2sce(ts_largeintestine_h5ad)
#
sce.ts_intestine
table(sce.ts_intestine$cell_type)
table(sce.ts_intestine$broad_cell_class)
table(sce.ts_intestine$donor_id, sce.ts_intestine$broad_cell_class )

rowData(sce.ts_intestine)


# Not needed, but first filter down to matched genes in our panel
sce.ts_intestine.genename <- sce.ts_intestine[rowData(sce.ts_intestine)$feature_name %in% rownames(so),]

# are there any duplicates (we'd need to handle them)
sum(table(rowData(sce.ts_intestine.genename)$feature_name) != 1 )

# just rename the genes to the gene names
rownames(sce.ts_intestine.genename) <-  rowData(sce.ts_intestine.genename)$feature_name



# Pull out the normalised matrix.
# Quirk of this coming from the python world, there's no counts assay, and the normalised is in 'X'
ref_norm_matrix <- assay(sce.ts_intestine.genename, 'X')



# Show a heatmap of htose at the sample summary level?

norm_matrix <- GetAssayData(so, assay = 'RNA', layer = 'data')

predictions <- SingleR::SingleR(test = norm_matrix,
                                ref   = ref_norm_matrix,
                                labels = sce.ts_intestine.genename$broad_cell_class)




so$singleR_labels <- predictions$labels
DimPlot(so, group.by='singleR_labels')
