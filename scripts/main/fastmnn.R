
working_dir <- '/data/ShernData/NF1TumorAtlas/'; setwd(working_dir)

source('scripts/main/load_packages.r')

object.list <- readRDS(snakemake@input[['seurat_object']])


object.list <- sapply(object.list, CleanVarGenes, nHVG=3000, simplify=FALSE)

features <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

object <- object.list %>% RunFastMNN(k=20, features=features, auto.merge=TRUE)

Project(object) <- 'NF1_mnn_integrate'


saveRDS(object, snakemake@output[['seurat_object']])

