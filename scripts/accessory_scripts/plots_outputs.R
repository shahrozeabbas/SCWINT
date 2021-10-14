plots_outputs <- function(object, out_dir=getwd(), w=20, h=12, group='seurat_clusters', plot=FALSE, output=TRUE, markers=TRUE) {
  
  # plot <- FALSE; output <- FALSE
  require(data.table); require(Seurat); require(ggplot2); require(dplyr); require(SeuratDisk)

  to_disk <- '/gpfs/gsfs12/users/ShernData/'
  root <- paste0(out_dir, '/final_', Project(object))

  if (markers) {
    system(paste('sbatch', paste0('--job-name=', Project(object), '_markers'), paste0(to_disk, '/job_submitters/markers.sh'), root, group))
  }
  
  
  if (plot) output <- TRUE

  if (output) {

    md <- copy(object@meta.data)
    setDT(md, keep.rownames=TRUE)
    setnames(md, old=c('rn', 'orig.ident', 'nCount_RNA', 'nFeature_RNA'), 
            new=c('cells', 'scafID', 'nUMI', 'nGene'))
    
    umap <- copy(as.data.frame(object[['umap']]@cell.embeddings))
    setDT(umap, keep.rownames='cells')
    
    md <- md[umap, on='cells']; rm(umap)
    fwrite(md, paste0(root, '_metadata.csv'))


    avg <- as.data.frame(object %>% AverageExpression(group.by=group, assay=DefaultAssay(object)))
    setDT(avg, keep.rownames='genes')

    fwrite(x=avg, file=paste0(root, '_cluster_average_expression.csv'))

    genes <- avg[, genes]
    avg[, genes := NULL]
    setDF(avg, rownames=genes)
    
    # adj <- WGCNA::adjacency(datExpr=avg, type='signed', power=1)
    # hm_colors <- colorRampPalette(wesanderson::wes_palette('Zissou1', type='discrete'))
    # pheatmap::pheatmap(mat=adj, show_rownames=TRUE, show_colnames=TRUE, 
    #         filename=paste0(root, '_cluster_average_correlation.pdf'), 
    #         color=hm_colors(100))


    fwrite(x=md[, {
                totwt=.N
                .SD[, .(percent=.N / totwt), by=pathology]
                      }, by=seurat_clusters][md[, .(nCells=.N), keyby=.(seurat_clusters, pathology)], 
                                      on=c('seurat_clusters', 'pathology')], 
          file=paste0(root, '_cluster_pathology_composition.csv')) 
    
    
    if (plot) {

      g <- ggplot(data=md, aes(x=UMAP_1, y=UMAP_2)) + theme_classic()

      pal <- wesanderson::wes_palette('Zissou1', type='discrete')[c(1, 3, 5)]
      continuous_stats <- c('nUMI', 'nGene', 'percent.mt', 'doublet_scores')
      
      invisible(parallel::mclapply(X=continuous_stats, FUN=function(attribute) {
        
        middle <- md[, mean(c(max(get(attribute)), min(get(attribute))))]
        ggsave(plot=g +
                geom_point(aes(color=get(attribute)), alpha=0.7, size=0.01) + 
                scale_color_gradient2(low=pal[1], mid=pal[2], high=pal[3], midpoint=middle) +
                labs(color=attribute),
              filename=paste0(out_dir, '/', attribute, '_umap.pdf'), width=w, height=h)
        
      }, mc.cores=length(continuous_stats)))
        
      
      pal <- colorRampPalette(ggsci::pal_ucscgb()(26))
      discrete_stats <- c('pathology', 'predicted_doublets', 'Phase', 'seurat_clusters', 'specimen', 'project')
      
      invisible(parallel::mclapply(X=discrete_stats, FUN=function(attribute) {

        n_colors <- nrow(md[, .N, keyby=get(attribute)])  
        ggsave(plot=g + 
                geom_point(aes(color=get(attribute)), size=0.01) + 
                scale_color_manual(values=pal(n_colors)) +
                labs(color=attribute) + 
                guides(color=guide_legend(override.aes=list(size=5))), 
              filename=paste0(out_dir, '/', attribute, '_umap.pdf'), width=w, height=h)
        
      }, mc.cores=length(discrete_stats)))
    }
  }
  
}


