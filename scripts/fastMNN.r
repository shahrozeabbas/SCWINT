args <- commandArgs(trailingOnly=TRUE)

pkgs <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'future', 'SeuratWrappers', 'SeuratDisk')
invisible(sapply(pkgs, require, character.only=TRUE))

home <- '/insert/main/working/directory/'
data_dir <- paste0(home, 'data2integrate/')
if (!dir.exists(data_dir)) stop('data directory does not exist!')

setwd(data_dir)
output_dir <- paste0(home, 'fastMNN_output/')
if (!dir.exists(output_dir)) dir.create(output_dir)

project_list_name <- 'final_NF1_project_list.rds'
project_list <- readRDS(project_list_name)

accessory_scripts <- list.files(path=paste0(home, 'accessory_scripts'), pattern='.R', full.names=TRUE)
invisible(lapply(X=accessory_scripts, FUN=source))


k <- as.numeric(args[1])
roundUP <- function(x, m) {x + m - x %% m}
plan('multicore', workers=roundUP(length(project_list), 2))
options(future.globals.maxSize=100 * 1000 * 1024^2)


project_list <- lapply(project_list, function(object) {
    object %>%
        filter_cells() %>%
        NormalizeData() %>%
        FindVariableFeatures()
})

if (k > 1) {
    combined <- project_list %>% RunFastMNN(k=k, auto.merge=TRUE)
    neighbors <- k / 2
}

if (k > 0 & k < 1) {
    combined <- project_list %>% RunFastMNN(prop.k=k, auto.merge=TRUE)
}

output_dir <- ifelse(k < 1,
                paste0(output_dir, 'anchor', gsub('.', '', k, fixed=TRUE)),
                    paste0(output_dir, 'anchor', k))

if (!dir.exists(output_dir)) dir.create(output_dir)

combined <- combined %>%
                    RunUMAP(reduction='mnn', dims=seq(30)) %>%
                    FindNeighbors(reduction='mnn', dims=seq(30)) %>%
                    FindClusters(resolution=1.5) %>%
                    annotate_clusters(path2plots=output_dir, assay='RNA')

Project(combined) <- 'NF1_mnn_integrate'

saveRDS(combined, paste0(output_dir, '/NF1_fastMNN_integrate.rds'))
combined %>% SaveH5Seurat(filename=paste0(output_dir, '/final_', Project(combined), '.h5Seurat'), overwrite=TRUE)

# DefaultAssay(combined) <- 'mnn.reconstructed'
combined %>% plots_outputs(out_dir=output_dir, markers=TRUE, plot=TRUE)
