
PKGS <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'SeuratWrappers', 'parallel', 'future', 'reticulate', 'SoupX', 'scales')

invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files('scripts/utility', full.names=TRUE, pattern='\\.r$'), source))


data_path <- 'data/'

setDTthreads(threads=1)
