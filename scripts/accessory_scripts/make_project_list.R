make_project_list <- function(project_path=getwd(), projects=NULL, plot_doublets=TRUE, verbose=TRUE, output_dir=getwd()) {
  
  # n_projects <- length(list.dirs(recursive=FALSE))
  # ncores <- ifelse(availableCores() > n_projects, n_projects, availableCores() / 2)
  # read in master table
  if (is.null(projects)) {
    projects <- fread(paste0(project_path, '/NF1_project_scaf_cellcount.csv'))
  } else {
    projects <- fread(projects)
  }
  # begin looping through datasets
  project_list <- sapply(X=list.dirs(recursive=FALSE, full.names=FALSE), FUN=function(project) {

    if (verbose) print(project)
    
    setwd(paste0(project_path, '/', project, '/SCAFS/'))
  
    if (length(list.dirs(recursive=FALSE)) != 0) {
    
      scaf_list <- sapply(X=list.dirs(recursive=FALSE, full.names=FALSE), FUN=function(scaf) {
        
        if (verbose) print(scaf)
        
        data_dir <- paste0(scaf, '/outs/filtered_feature_bc_matrix.h5')
        seurat <- CreateSeuratObject(counts=Read10X_h5(filename=data_dir),
                                    project=scaf, min.cells=0, min.features=0)
        
        # add mitochondrial percentages
        seurat[['percent.mt']] <- PercentageFeatureSet(seurat, pattern='^MT-')
        
        # grab pathology based on SCAF for metadata
        pathology <- projects[ProjectID %chin% project & ScafID %chin% tstrsplit(scaf, '_')[[1]], 
                      rep(Pathology, nCells)]
        specimen <- projects[ProjectID %chin% project & ScafID %chin% tstrsplit(scaf, '_')[[1]], 
                      rep(SpecimenID, nCells)]

        
        names(specimen) <- names(pathology) <- seurat$orig.ident
        # add pathology as meta data
        seurat <- seurat %>% 
                    AddMetaData(metadata=factor(pathology), col.name='pathology') %>%
                    AddMetaData(metadata=factor(specimen), col.name='specimen')

        
        doublet_rate <- (ncol(seurat) / 1000) * 0.008
        
        seurat <- seurat %>% 
                              scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate, exp_rate_call=TRUE,
                                      plot=ifelse(plot_doublets, paste0(output_dir, '/', scaf, '_'), FALSE), 
                                      remove_doublets=TRUE)
                              
                              # doublet_finder(npcs=30, expected_doublet_rate=doublet_rate, cores=ncores) %>% 

                              # scds(expected_doublet_rate=doublet_rate)                                        
         
                  
          })
        
        }

      if (length(scaf_list) == 1) {
        merge_scafs <- scaf_list[[1]]
        Project(merge_scafs) <- project
      
      } else {
        merge_scafs <- merge(x=scaf_list[[1]], y=scaf_list[-1], 
                                    add.cell.ids=tstrsplit(names(scaf_list), '_')[[1]], project=project)
        }

      p <- rep(project, ncol(merge_scafs))
      names(p) <- merge_scafs$orig.ident
      merge_scafs <- merge_scafs %>% AddMetaData(metadata=factor(p), col.name='project')  
    
    })

  return(project_list[which(sapply(project_list, class) == 'Seurat')])

}



split_project_list <- function(project_list, split=NULL) {
  
  project_list <- unlist(sapply(X=project_list, FUN=function(project) {
    project %>% SplitObject(split.by=split)
  })) 

  names(project_list) <- tstrsplit(names(project_list), '.', fixed=TRUE)[[2]]

  for (counter in seq_along(project_list)) {
    Project(project_list[[counter]]) <- names(project_list)[counter]
  }
  
  return(project_list)
}
