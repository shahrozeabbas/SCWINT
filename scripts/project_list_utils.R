args <- commandArgs(trailingOnly=TRUE)

pkgs <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'reticulate', 'SeuratDisk', 'future')
invisible(sapply(pkgs, require, character.only=TRUE))

root <- '/gpfs/gsfs12/users/ShernData/'

use_python(system('which python', intern=TRUE), required=TRUE)
source_python(paste0(root, 'scripts_NF1/scrublet_py.py'))


name <- 'NF1_merged'
merge2all <- FALSE

source(paste0(root, 'scripts_NF1/scds.R'))
source(paste0(root, 'scripts_NF1/scrublet.R'))
source(paste0(root, 'scripts_NF1/doublet_finder.R'))
source(paste0(root, 'scripts_NF1/make_project_list.R'))

setwd(paste0(root, 'NF1TumorData/'))
project_list <- make_project_list(plot_doublets=FALSE)

if (merge2all) {
    merge_projects <- merge(x=project_list[[1]], y=project_list[-1], 
                           add.cell.ids=names(project_list), project=name)
    
    merge_projects %>% SaveH5Seurat(filename=paste0(root, 'NF1TumorData/NF1_merged.h5Seurat', overwrite=TRUE))                    
} 


project_list <- project_list %>% split_project_list(split='specimen')
saveRDS(project_list, file=paste0(root, 'NF1TumorData/', args[1], '.rds'))


# pnf <- lapply(project_list, function(project) if (unique(project$pathology) == 'PN') project)
# anf <- lapply(project_list, function(project) if (unique(project$pathology) == 'ANNUBP') project)
# mpnst <- lapply(project_list, function(project) if (unique(project$pathology) == 'MPNST') project)


# pnf <- pnf[lengths(pnf) != 0]
# anf <- anf[lengths(anf) != 0]
# mpnst <- mpnst[lengths(mpnst) != 0]

# types <- c('pnf', 'anf', 'mpnst')

# for (type in types) {

#     saveRDS(eval(as.symbol(type)), 
#             file=paste0(root, 'NF1TumorData/', type, '_project_list.rds'))
    
# }
