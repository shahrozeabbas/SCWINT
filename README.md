# A Single Cell RNA-Seq Workflow to Integrate and Normalize Diverse Datasets in NF1 Nerve Tumors


## Background
Single cell sequencing allows researchers to explore disease mechanisms with greater resolution than microarray or bulk-tissue. One drawback is that data is sparse and difficult to interpret. Recent approaches have leveraged the integration of multiple datasets and modalities to further elucidate cell types and discover populations that would otherwise be difficult to observe. Several methods exist for data integration. Here, we combine expression information of 49 different tumors from different reagents and chemistries using FastMNN to normalize and smooth the data for visualization and differential gene expression analysis.  
![image](https://user-images.githubusercontent.com/28969387/168441675-a25aad12-a771-423d-8a8f-ddc3e6c7d7c9.png)

## Workflow and Methods
FastMNN’s uses a method of matching mutual nearest neighbors to remove batch effect. Matching MNN pairs requires the user to properly select the number of neighbors to be matched. Samples were sequenced to create count matrices representing expression information of every gene across ~400,000 cells. Ambient RNA was removed using SoupX. Doublets were identified using Scrublet and removed using an expected rate from 10X Genomics. Mutual nearest neighbors were matched using 20 cells per sample. Louvain community detection for clustering and marker identification was performed using Seurat’s native methods. In the flowchart below, preprocessing and QC steps are labeled in yellow, the main integration steps are in blue, and annotation and subtype analysis in red.
![image](https://user-images.githubusercontent.com/28969387/168441695-f594ffbe-e9e5-4251-adcf-3f0aee8127cf.png)

<img width="775" alt="image" src="https://user-images.githubusercontent.com/28969387/168441700-b3245e01-c513-4d71-b0e1-449b13359df1.png">

## Ambient RNA Removal with SoupX
This figure from the SoupX paper outlines the process of correcting count matrices for cell-free RNA captured by droplet based techniques. This additional ‘soup’ adds to transcriptomic profile of a cell, negatively impacting downstream analyses. 
![image](https://user-images.githubusercontent.com/28969387/168441812-ad20c9af-383e-4a80-a231-efffef057dc2.png)
<img width="766" alt="image" src="https://user-images.githubusercontent.com/28969387/168441760-e1c9fa83-f381-42a3-b0fb-c4cbf232e5a3.png">

## Doublet Removal with Scrublet
Doublet scores calculated by Scrublet. Green bars indicate scores for singlet cells, while purple bars indicate potential doublets. Presence of doublets has been shown to confound cell type clustering.
![image](https://user-images.githubusercontent.com/28969387/168441867-b5bd629c-f7b3-4234-9e17-24bfc7ece963.png)
<img width="751" alt="image" src="https://user-images.githubusercontent.com/28969387/168441870-8358e11e-ae21-43cf-a526-a8583d925d3e.png">







A `yml` file is included to recreate the [conda](https://www.anaconda.com) environment used for the analysis by running:   
`conda env create -f environment.yml`
