working_dir <- '/data/ShernData/NF1TumorAtlas/'; setwd(working_dir)

source('scripts/main/load_packages.r')

# set soupx contamination rate
contamination_rate <- 0.20
reticulate::source_python('scripts/utility/scrublet_py.py')

# read in master table
projects <- fread(snakemake@input[['samples']])

# begin looping through datasets
project_list <- sapply(X=projects[, .N, keyby=ProjectID][[1]], FUN=function(project) {
  
  project_path <- paste0(data_path, project, '/')

  scaf_list <- sapply(X=list.dirs(project_path, recursive=FALSE, full.names=FALSE), FUN=function(scaf) {

    if (scaf %chin% projects[, ScafID]) {

      lane_dir <- paste0(project_path, scaf, '/outs/')

      suppressWarnings({
        message('creating SoupX corrected seurat object...')

        sc <- load10X(lane_dir)
        sc <- setContaminationFraction(sc, contamination_rate)
        out <- adjustCounts(sc)
        
        seurat <- CreateSeuratObject(counts=out, project=scaf, min.cells=0, min.features=0)

      })

      # add mitochondrial percentages
      seurat[['percent.mt']] <- PercentageFeatureSet(seurat, pattern='^MT-')
      seurat[['percent.rb']] <- PercentageFeatureSet(seurat, pattern='^RP[SL]')

      m <- copy(seurat@meta.data)
      setDT(m, keep.rownames='cells')

      # grab pathology based on LANE for metadata
      m[, `:=` (
        pathology=projects[ProjectID %chin% project & ScafID %chin% scaf, Pathology],
        specimen=projects[ProjectID %chin% project & ScafID %chin% scaf, SpecimenID],
        tumor=projects[ProjectID %chin% project & ScafID %chin% scaf, TumorID],
        scaf=projects[ProjectID %chin% project & ScafID %chin% scaf, ScafID]
      )]

      scaf <- m[, scaf]; tumor <- m[, tumor]; specimen <- m[, specimen]; pathology <- m[, pathology]

      names(scaf) <- names(specimen) <- names(pathology) <- names(tumor) <- m[, cells]
      # add pathology as meta data

      seurat <- seurat %>% 
                  AddMetaData(metadata=factor(pathology), col.name='pathology') %>%
                  AddMetaData(metadata=factor(specimen), col.name='specimen') %>% 
                  AddMetaData(metadata=factor(tumor), col.name='tumorID') %>% 
                  AddMetaData(metadata=factor(scaf), col.name='laneID')

      
      doublet_rate <- (ncol(seurat) / 1000) * 0.008

      seurat <- seurat %>% scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate)

      }

  })
  
  scaf_list <- scaf_list[which(sapply(scaf_list, class) == 'Seurat')]


  if (length(scaf_list) == 1) {

    merge_scafs <- scaf_list[[1]]
    Project(merge_scafs) <- project
  
  } else {
    merge_scafs <- merge(x=scaf_list[[1]], y=scaf_list[-1], add.cell.ids=names(scaf_list), project=project)
  }

  p <- rep(project, ncol(merge_scafs))
  names(p) <- merge_scafs$orig.ident
  merge_scafs <- merge_scafs %>% AddMetaData(metadata=factor(p), col.name='project')  

})

project_list[which(sapply(project_list, class) == 'Seurat')]


object <- merge(x=project_list[[1]], y=project_list[-1], add.cell.ids=names(project_list), project='tumor_atlas')

metadata <- copy(object@meta.data)
setDT(metadata, keep.rownames='cells')
fwrite(metadata, file=snakemake@output[['metadata']])

saveRDS(object, snakemake@output[['seurat_object']])

