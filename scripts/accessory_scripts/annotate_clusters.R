
    
annotate_clusters <- function(object, path2plots=getwd(), h=12, w=20, assay=DefaultAssay(object)) {

    fibroblast <- c('DCN', 'APOD', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP2', 'TAGLN', 'THY1')

    CTL <- c('CD3D', 'CD3E', 'CD8A', 'CD8B', 'GZMB', 'CCL5')

    monocyte <- c('LYZ', 'S100A9', 'S100A8')

    macrophage <- c('FOLR2', 'CD14', 'C1QA', 'C1QB', 'C1QC')

    endothelial <- c('VWF', 'PLVAP', 'HSPG2', 'MMRN1', 'LYVE1')

    SC <- c('S100B', 'CAPS', 'MIA', 'CDH19', 'LGI4', 'PLP1')

    mast <- c('TPSAB1', 'HPGDS', 'CTSG', 'KIT')

    pericyte <- c('RGS5', 'PDGFRB', 'CCDC102B', 'NOTCH3', 'TINAGL1')

    MPNST <- c('DLK1', 'CYTL1', 'SNAI2', 'IGF2')

    object <- object %>% 
                AddModuleScore(features=fibroblast, ctrl=100, name='Fibroblast', assay=assay) %>%
                AddModuleScore(features=CTL, ctrl=100, name='CTL', assay=assay) %>%
                AddModuleScore(features=monocyte, ctrl=100, name='Monocyte', assay=assay) %>% 
                AddModuleScore(features=macrophage, ctrl=100, name='Macrophage', assay=assay) %>% 
                AddModuleScore(features=endothelial, ctrl=100, name='Endothelial', assay=assay) %>%
                AddModuleScore(features=SC, ctrl=100, name='SC', assay=assay) %>%
                AddModuleScore(features=pericyte, ctrl=100, name='Pericyte', assay=assay) %>%
                AddModuleScore(features=MPNST, ctrl=100, name='MPNST', assay=assay)

    p <- object %>% 
        FeaturePlot(features=c('Fibroblast1', 'CTL1', 'Monocyte1', 'Endothelial1', 'SC1', 'Macrophage1', 'Pericyte1', 'MPNST1'), 
        min.cutoff='q10', max.cutoff='q90', ncol=3, cols=c('lightgrey', 'red'))
    
    ggsave(plot=p, height=h, width=w, filename=paste0(path2plots, paste0('/NF1_marker_umap.pdf')))
    
    return(object)
    
}
    
    
