args <- commandArgs(trailingOnly=TRUE)
type <- args[1]; action <- args[2]

message(type); message(action)

pkgs <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'future', 'SeuratDisk')
invisible(sapply(pkgs, require, character.only=TRUE))

to_disk <- '/gpfs/gsfs12/users/ShernData/'
output_dir <- paste0(to_disk, type, '_integrate')
output_dir <- paste0(output_dir, '/default_scrubbed_specimen_merge')
if (!dir.exists(output_dir)) dir.create(output_dir)


setwd(paste0(to_disk, 'NF1TumorData/'))

invisible(lapply(X=list.files(path=paste0(to_disk, 'scripts_NF1'),
                                pattern='.R', full.names=TRUE), FUN=source))


setDTthreads(threads=2)
plan('multicore', workers=availableCores() / 2)


filter <- cluster <- FALSE
if (action == 'filter') filter <- TRUE
if (action == 'cluster') cluster <- TRUE
if (action == 'both') filter <- cluster <- TRUE

d <- 50
hvg <- 3000
ngbs <- 200
groupby <- 'orig.ident'
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


project_list_name <- 'NF1_project_list_scrubbed_specimen_merge.rds'
message(project_list_name)


if (filter) {

       project_list <- readRDS(file=project_list_name)

       project_list <- sapply(X=project_list, FUN=function(project) {

              project <- project %>%
                     qc(nhvg=hvg, group=groupby, score_cc=FALSE, elbow=FALSE, dims=d,
                            path2plots=paste0(output_dir, '/unfiltered/unfiltered_', Project(project))) %>%

                     filter_cells() %>% qc(nhvg=hvg, group=groupby, score_cc=FALSE, elbow=FALSE, dims=d,
                                          path2plots=paste0(output_dir, '/filtered/filtered_', Project(project)), write_plots=TRUE)

       })


       features <- SelectIntegrationFeatures(object.list=project_list, nfeatures=hvg)
       save(features, file=paste0(output_dir, '/integration_features_', hvg, '_', type, '.RData'))

       project_list <- sapply(X=project_list, FUN=function(project) {
              project <- project %>%
                     qc(group=groupby, features=features, regress_noise=TRUE, normalize=FALSE, elbow=FALSE,
                     dims=d, path2plots=paste0(output_dir, '/filtered/filtered_regressed_', Project(project)), write_plots=TRUE)
       })

       saveRDS(object=project_list, file=paste0(output_dir, '/filtered_', project_list_name))

}

if (cluster) {

       k <- as.numeric(args[3])
       # r <- as.numeric(args[4])

       load(file=paste0(output_dir, '/integration_features_', hvg, '_', type, '.RData'))
       project_list <- readRDS(file=paste0(output_dir, '/filtered_', project_list_name))

       anchors <- FindIntegrationAnchors(object.list=project_list, anchor.features=features, reduction='rpca', k.anchor=k, scale=FALSE)

       # this command creates an 'integrated' data assay
       combined <- IntegrateData(anchorset=anchors)

       DefaultAssay(combined) <- 'integrated'

       output_dir <- paste0(output_dir, '/anchor', k)
       if (!dir.exists(output_dir)) dir.create(output_dir)

       # Run the standard workflow for visualization and clustering
       combined <- combined %>% qc(group='specimen', score_cc=FALSE, features=features, dims=30,
                                   normalize=FALSE, path2plots=paste0(output_dir, '/filtered/integrated_'),
                                   write_plots=TRUE)


       combined <- combined %>%
              FindNeighbors(reduction='pca', dims=seq(30), k.param=100) %>% FindClusters() %>%
              RunUMAP(reduction='pca', dims=seq(30)) %>% annotate_clusters(path2plots=output_dir, assay='RNA')

       Project(combined) <- paste0(type, '_integrated')

       combined %>% SaveH5Seurat(filename=paste0(output_dir, '/final_', Project(combined), '.h5Seurat'), overwrite=TRUE)

       saveRDS(combined, paste0(output_dir, '/final_', Project(combined), '.rds'))

       # for (resolution in 6:12 / 10) system(paste('sbatch add_resolutions.sh', resolution))
       combined %>% plots_outputs(out_dir=output_dir, markers=TRUE, plot=TRUE)
       message(paste('Finished:', type, action, k))

}
