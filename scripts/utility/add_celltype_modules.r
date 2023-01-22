
    
AnnotateClusters <- function(object, plot.out=getwd(), h=9, w=12, assay=DefaultAssay(object), cell.size=0.1) {

    genes <- list(
        Fibroblast = c('DCN', 'APOD', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP2', 'TAGLN', 'THY1'),
        CTL = c('CD3D', 'CD3E', 'CD8A', 'CD8B', 'GZMB', 'CCL5'),
        Monocyte = c('LYZ', 'S100A9', 'S100A8'),
        Macrophage = c('FOLR2', 'CD14', 'C1QA', 'C1QB', 'C1QC'),
        Endothelial = c('VWF', 'PLVAP', 'HSPG2', 'MMRN1', 'LYVE1'),
        SC = c('S100B', 'CAPS', 'MIA', 'CDH19', 'LGI4', 'PLP1'),
        Mast = c('TPSAB1', 'HPGDS', 'CTSG', 'KIT'),
        Pericyte = c('RGS5', 'PDGFRB', 'CCDC102B', 'NOTCH3', 'TINAGL1'),
        MPNST = c('DLK1', 'CYTL1', 'SNAI2', 'IGF2')
    )

    colors <- c('lightgrey', 'navy')
    object <- object %>% AddModuleScore(features=genes, assay=assay, name=names(genes))

    p <- object %>% 
        FeaturePlot(features=paste0(names(genes), seq(length(genes))), pt.size=cell.size, 
        min.cutoff='q10', max.cutoff='q90', ncol=3, cols=colors, reduction='umap')
    
    ggsave(plot=p, height=h, width=w, filename=paste0(plot.out, 'activity_module_umap.pdf'))
    
    return(object)
    
}
    
    