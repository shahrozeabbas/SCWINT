cell_cycle_qc <- function(object = NULL, write_plots = FALSE, path2plots = NULL) {
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  noise <- c('S.Score', 'G2M.Score')
  cc_plots <- vector(mode = 'list')
  
    
  object <- object %>% 
    ScaleData(features = rownames(object)) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
 
  object <- object %>% RunPCA(features = c(s.genes, g2m.genes))
  
  cc_pca_before <- object %>% DimPlot()
  
  object <- object %>%
    ScaleData(vars.to.regress = noise, features = rownames(object))
  
  object <- object %>% RunPCA(features = c(s.genes, g2m.genes))
  
  cc_pca_after <- object %>% DimPlot()
  
  cc_plots[['cc_before']] <- cc_pca_before
  cc_plots[['cc_after']] <- cc_pca_after

  if (write_plots) {
    
    invisible(lapply(X = seq_along(cc_plots), FUN = function(counter) {
      ggsave(plot = cc_plots[[counter]], width = 10, height = 8, 
             filename = paste0(path2plots, '_qc_plot_', names(cc_plots)[counter], '.pdf'))
    }))
  }
  
  return(object)
}


