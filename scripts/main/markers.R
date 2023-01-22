
working_dir <- '/data/ShernData/NF1TumorAtlas/'; setwd(working_dir)

source('scripts/main/load_packages.r')

ngbs <- 50
plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

object <- readRDS(snakemake@input[['seurat_object']])

Idents(object) <- 'seurat_clusters'

g <- data.table(genes=rownames(object[['RNA']]))
clean.genes <- g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]

markers <- object %>% FindAllMarkers(assay='RNA', min.pct=0.25, features=clean.genes)

fwrite(x=markers, file=snakemake@output[['markers']])


metadata <- copy(object@meta.data)
setDT(metadata, keep.rownames='cells')

umap <- copy(as.data.frame(object[['umap']]@cell.embeddings))
setDT(umap, keep.rownames='cells')

metadata <- metadata[umap, on='cells']
fwrite(metadata, file=snakemake@output[['metadata']])


groups <- c('seurat_clusters', 'tumorID', 'pathology', 'project')
titles <- c('Clusters', 'Tumors', 'Pathology', 'Project')

plot.list <- lapply(seq(4), function(counter) {
    object %>% 
        DimPlot(group.by=groups[counter], label=TRUE, repel=TRUE, reduction='umap', pt.size=0.05) + 
            ggtitle(titles[counter]) + xlab('') + ylab('') + NoLegend()
})

g <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=2, nrow=2)

ggsave(plot=g, width=18, height=11, filename=snakemake@output[['umap']])
