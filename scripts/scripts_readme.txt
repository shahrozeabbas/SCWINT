
`fastMNN.r` integrates samples using the mutual nearest neighbors approach.

`integration_NF1.r` integrates samples using Seurat's rPCA method.

`project_list_utils.R` is used to create a list to integrate after merging between experiments. 

`run_markers.r` reads an h5 Seurat object file then finds markers per cluster using Seurat's `FindAllMarkers` at default. 
