qc <- function(object=NULL, group=NULL, split=NULL, features=NULL, normalize=TRUE,
               regress_noise=FALSE, nhvg=2000, write_plots=TRUE, score_cc=TRUE, elbow=TRUE,
               jackstraw=FALSE, nPerm=100, w=14, h=10, path2plots=NULL, dims=100) {
  
  require(Seurat); require(data.table)
  
  noise <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
  
  if (class(object) == 'Seurat') {
    
    if (normalize) {
      object <- object %>% 
        NormalizeData() %>% 
        FindVariableFeatures(nfeatures=nhvg)
    }
    
    if (score_cc) {
      noise <- c(noise, 'S.Score', 'G2M.Score')
      object <- object %>% 
        cell_cycle_qc(write_plots=write_plots, path2plots=path2plots)
    }

    if (is.null(features)) {
      features <- VariableFeatures(object)
    } 
    
    if (regress_noise) {  
      object <- object %>% 
        ScaleData(vars.to.regress=noise, features=features, split.by=split)
    } else {
      object <- object %>% ScaleData(features=features, split.by=split) 
    }

    object <- object %>% RunPCA(features=features, npcs=dims)
    
    if (write_plots) {
      q <- object %>% VlnPlot(features=noise[1:3], ncol=3, group.by=group, split.by=split)
      p <- object %>% DimPlot(group.by=group, split.by=split, reduction='pca')
      
      plot1 <- object %>% FeatureScatter(feature1=noise[2], feature2=noise[3], group.by=group)
      plot2 <- object %>% FeatureScatter(feature1=noise[2], feature2=noise[1], group.by=group)
      
      plot_list <- list(plot1 + plot2, q, p)
      
      if (jackstraw) {
        object <- object %>% JackStraw(num.replicate=nPerm) %>% ScoreJackStraw(dims=1:20)
        object %>% JackStrawPlot(dims=1:20) %>% ggsave(width=w, height=h, 
                                                    filename=paste0(path2plots, '_js_plot.pdf'))
        
      }

      if (elbow) {
        object %>% 
          ElbowPlot(ndims=dims) %>% 
          ggsave(width=w, heigh=h, filename=paste0(path2plots, '_elbow_plot.pdf'))
      }
      
  
      invisible(lapply(X=seq_along(plot_list), FUN=function(counter) {
        ggsave(plot=plot_list[[counter]], width=w, height=h, 
               filename=paste0(path2plots, '_qc_plot_', counter, '.pdf'))
             }))
    }
    
  }
  
  return(object)
}



