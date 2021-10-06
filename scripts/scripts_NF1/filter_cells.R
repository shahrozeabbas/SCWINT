filter_cells <- function(object, lower=0.02, upper=0.98) {
  require(data.table); require(Seurat)

  metadata <- copy(object@meta.data)
  
  setDT(metadata, keep.rownames=TRUE)
  new_names <- c('cells', 'scafID', 'nUMI', 'nGene')
  setnames(metadata, old=colnames(metadata)[1:4], new=new_names)
  
  cells2keep <- metadata[fintersect(fintersect(metadata[, 
                                    .(cellNo=.I[between(nUMI, quantile(nUMI, lower), quantile(nUMI, upper))]), 
                                    by=scafID], 
                           metadata[, 
                                    .(cellNo=.I[between(nGene, quantile(nGene, lower), quantile(nGene, upper))]), 
                                    by=scafID]), metadata[, 
                                                            .(cellNo=.I[between(percent.mt, quantile(percent.mt, lower), quantile(percent.mt, upper))]), 
                                                            by=scafID])[, cellNo]][nGene > 300 & nGene < 6000, cells]
  
  
  
  return(object %>% 
    subset(cells=cells2keep))

}



