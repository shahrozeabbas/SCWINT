working_dir <- '/data/ShernData/NF1TumorAtlas/'; setwd(working_dir)

source('scripts/main/load_packages.r')

object <- readRDS(snakemake@input[['seurat_object']])

DefaultAssay(object) <- 'mnn.reconstructed'

object <- object %>% 
    RunUMAP(dims=1:30, reduction='mnn') %>%
    FindNeighbors(dims=1:30, reduction='mnn') %>% 
    FindClusters(resolution=0.8, algorithm=1)

saveRDS(object, snakemake@output[['seurat_object']])