scrublet <- function(object, min_counts=2, min_cells=3, expected_doublet_rate=0.06, min_gene_variability_pctl=85, 
                     n_prin_comps=50, sim_doublet_ratio=2, n_neighbors=NULL, source_scrublet=FALSE, plot=FALSE, 
                     remove_doublets=TRUE, exp_rate_call=TRUE) 
{
  if (!py_available('scrublet')) stop('python module scrublet does not seem to be installed; - try running "py_config()"')
  
  if (class(object) != 'Seurat') stop('object is not of type "Seurat"') 
  
  if (source_scrublet) {
    scrublet_path <- '/gpfs/gsfs12/users/ShernData/scripts_NF1/scrublet_py.py'
    source_python(scrublet_path)
  }
  
  X <- as(t(as.matrix(GetAssayData(object=object, slot='counts'))), 'TsparseMatrix')
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
  dim <- as.integer(X@Dim)
  
  if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(nrow(X)))
  
  
  scrublet_py_args <- c(list(i=i, j=j, val=val, dim=dim,
                             expected_doublet_rate=expected_doublet_rate, min_counts=min_counts, 
                             min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps, 
                             sim_doublet_ratio=sim_doublet_ratio, n_neighbors=n_neighbors))
  
  scrublet_res <- do.call(scrublet_py, scrublet_py_args)
  names(scrublet_res) <- c('doublet_scores', 'predicted_doublets')
  
  d <- data.table(cells=object$orig.ident, doublet_score=scrublet_res$doublet_scores)
  
  doublet_scores <- d[, doublet_score]
  names(doublet_scores) <- d[, cells]
  
  if (exp_rate_call) {
    
    predicted_doublets <- d[, fifelse(
          doublet_score > quantile(doublet_score, 1 - expected_doublet_rate), 
                                                  'doublet', 'singlet')]

  } else {

    predicted_doublets <- d[, fifelse(
          doublet_score > max(mixtools::normalmixEM(doublet_score)$mu), 
                                                  'doublet', 'singlet')]
  }
  
  names(predicted_doublets) <- d[, cells]
  
  if (is.character(plot)) {
    p <- ggplot(data=d, aes(x=doublet_score, fill=cluster)) + 
      geom_density(alpha=0.5) + scale_fill_brewer(palette='Dark2') + theme_bw()

    p %>% ggsave(filename=paste0(plot, 'scrublet_scores.pdf'))
  }
  
  object <- object %>% AddMetaData(metadata=doublet_scores, col.name='doublet_scores') %>%
    AddMetaData(metadata=factor(predicted_doublets), col.name='predicted_doublets')

  if (remove_doublets) {
    m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
    singlets <- m[predicted_doublets %in% 'singlet', cells]
    object <- object %>% subset(cells=singlets)
  }

  return(object)
  
}

