

rule all:
    input:
        'output/filtered_final_metadata.csv',
        'output/seurat_scvi_integrated_cluster_markers.csv',
        'plots/seurat_mnn_seurat_clusters_tumors_pathology_project_umap.pdf'


rule preprocess:
    input: 
        samples='input/samples.csv'
    output:
        metadata='output/unfiltered_seurat_metadata.csv',
        seurat_object='objects/seurat_object_rna_preprocessed_01.rds'
    conda:
        'envs/singlecell.yml'
    script: 
        'scripts/main/preprocess.R'

rule filter_normalize:
    input:
        seurat_object='objects/seurat_object_rna_preprocessed_01.rds'
    output:
        seurat_object='objects/seurat_object_list_rna_preproccessed_filtered_normalized_02.rds'
    conda:
        'envs/singlecell.yml'
    threads:
        8
    script:
        'scripts/main/filter_normalize.R'

rule mnn:
    input:
        seurat_object='objects/seurat_object_list_rna_preproccessed_filtered_normalized_02.rds'
    output:
        seurat_object='objects/seurat_object_mnn_integrated_03.rds'
    conda:
        'envs/singlecell.yml'     
    script:
        'scripts/main/fastmnn.R'

rule umap_cluster:
    input:
        seurat_object='objects/seurat_object_mnn_integrated_03.rds'
    output:
        seurat_object='objects/seurat_object_mnn_integrated_clustered_04.rds'
    conda:
        'envs/singlecell.yml'
    script:
        'scripts/main/umap_cluster.R'

rule markers:
    input:
        seurat_object='objects/seurat_object_mnn_integrated_clustered_04.rds'
    output:
        metadata='output/filtered_final_metadata.csv',
        markers='output/seurat_mnn_integrated_cluster_markers.csv',
        umap='plots/seurat_mnn_seurat_clusters_tumors_pathology_project_umap.pdf'
    conda:
        'envs/singlecell.yml'
    threads:
        32
    script:
        'scripts/main/markers.R'

# rule annotate:
#     input:
#         seurat_object='objects/seurat_object_scvi_integrated_clustered_04.rds'
#     output:
#         annotation='output/sctype_seurat_clusters_annotation.csv'
#     conda:
#         'envs/singlecell.yml'
#     threads:
#         4
#     script:
#         'scripts/main/annotate_clusters.R'

# rule infercnv:
#     input:
#         reference='input/gencode_v19_gene_pos.txt',
#         annotation='output/sctype_seurat_clusters_annotation.csv',
#         seurat_object='objects/seurat_object_scvi_integrated_clustered_04.rds'
#     output:
#         seurat_object='objects/seurat_object_scvi_integrated_clustered_infercnv_05.rds'
#     params:
#         output_directory='output/infercnv'
#     conda:
#         'envs/singlecell.yml'
#     threads:
#         8
#     script:
#         'scripts/main/infercnv.R'