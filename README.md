# NF1-Integration

A `yml` file is included to recreate the [conda](https://www.anaconda.com) environment used for the analysis by running:   
`conda env create -f environment.yml`. 

Raw counts from each experiment were stored in [_Seurat_](https://satijalab.org/seurat/) objects for analysis. The list of objects to be integrated was created using `scripts/project_list_utils.R`. This script is used to add informative meta data such as doublet scores from [_Scrublet_](https://github.com/swolock/scrublet), and merge multiple sequencing lanes per experiment. Pipeline begins with calling `scripts/fastMNN.r`. Make sure paths are changed to reflect your workspace. 


The flowchart below outlines the integration pipeline using [_fastMNN_](https://marionilab.github.io/FurtherMNN2018/theory/description.html) from [SeuratWrappers](https://github.com/satijalab/seurat-wrappers). 


![fastMNN_1-3](https://user-images.githubusercontent.com/28969387/137391455-29614234-6615-4aa6-9c0a-d8ea438da26c.png)



