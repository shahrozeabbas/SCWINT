args <- commandArgs(trailingOnly=TRUE)

pkgs <- c('Seurat', 'dplyr', 'future')
invisible(sapply(pkgs, require, character.only=TRUE))


ngbs <- 10
plan('multicore', workers=availableCores() - 10)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

object <- SeuratDisk::LoadH5Seurat(file=paste0(args[1], '.h5Seurat'))

Idents(object) <- args[2]
markers <- object %>% FindAllMarkers(assay='RNA', min.pct=0.25)
data.table::fwrite(x=markers, file=paste0(args[1], '_cluster_markers.csv'))

