# miRTalk
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen?logo=github)](https://github.com/multitalk/miRTalk/actions) [![version>4](https://img.shields.io/badge/version-%3E%3D4.0-yellow?logo=Rstudio)](#install) 

### Extracellular vesicle-derived miRNA-mediated cell-cell communication inference for single-cell transcriptomic data

<img src='https://github.com/multitalk/miRTalk/blob/main/img/workflow.png'>

MicroRNAs are released from cells in extracellular vesicles (EVs), including exosomes and microvesicles, representing an essential mode of cell-cell communications via inhibitory effect on gene expression of receivers. The advent of single-cell RNA-sequencing (scRNA-seq) technologies has ushered in an era of elucidating EV-derived miRNA-mediated cell-cell communications at an unprecedented resolution. However, the lack of computational method to infer such communications from scRNA-seq data poses an outstanding challenge. Herein, we present miRTalk, a pioneering framework for discerning EV-derived miRNA-mediated cell-cell communications with a probabilistic model and a meticulously curated database, miRTalkDB, which catalogues EV-derived miRNA-target associations. Rigorous benchmarking against simulated and real-world datasets demonstrated the remarkable precision and robustness of miRTalk. Subsequently, we employed miRTalk to unravel the in-depth communicative mechanisms underlying three disease scenarios. In summary, miRTalk represents the first approach for inferring EV-derived miRNA-mediated cell-cell communications from scRNA-seq data, furnishing invaluable insights into the intercellular dynamics underpinning pathological processes.


# Install

```
# install devtools and install
install.packages(pkgs = "devtools"")
devtools::install_github("ZJUFanLab/miRTalk")
```

OR

```
# download the repository as ZIP
devtools::install_local("/path/to/miRTalk-main.zip")
```

# Usage
miRTalk method consists of three components, wherein the first is to infer the EV-derived miRNA across cells and the highly variable target genes, the second is to infer the cell-cell communication mediated by EV-derived miRNAs and their downstream targets. The third part is to visualize the miRNA-mediated cell-cell communication network and miRNA-target interaction network. Classification and description of miRTalk functions are shown in the document.

- ### Inference of EV-derived miRNA and highly variable target genes
```
# sc_data: A data.frame or matrix or dgCMatrix containing raw counts of single-cell RNA-seq data
# sc_celltype: A character containing the cell type of the single-cell RNA-seq data

> obj <- create_miRTalk(sc_data, sc_celltype, species = "Human", if_normalize = TRUE)
> obj
An object of class miRTalk 
0 EV-derived miRNA-target interactions

> obj <- find_miRNA(object = obj, mir_info = mir_info)
> obj <- find_hvtg(object = obj)
```

- ### Inference of cell-cell communication mediated by EV-derived miRNAs and their downstream targets
```
# object: miRTalk object after running find_miRNA() and find_hvtg()
# mir2tar: A data.frame containing the priori knowledge of miRNA-target interactions

> obj <- find_miRTalk(object = obj, mir2tar = mir2tar)
> obj
An object of class miRTalk 
2185 EV-derived miRNA-target interactions

> obj_cci <- get_miRTalk_cci(obj)
> str(obj_cci)
'data.frame':	2083 obs. of  9 variables:
 $ celltype_sender     : chr  "Bcell" "Bcell" "Bcell" "Bcell" ...
 $ celltype_receiver   : chr  "Bcell" "Bcell" "Bcell" "Bcell" ...
 $ miRNA               : chr  "hsa-miR-4426" "hsa-miR-29b-3p" "hsa-miR-29b-3p" "hsa-miR-29b-3p" ...
 $ miR_gene            : chr  "MIR4426" "MIR29B1" "MIR29B1" "MIR29B1" ...
 $ miRNA_activity      : num  0.429 0.481 0.481 0.481 0.481 ...
 $ target_gene         : chr  "PPIC" "CLDN1" "EREG" "TGFB3" ...
 $ target_gene_activity: num  5.81e-03 3.45e-04 1.03e-04 3.43e-04 4.85e-05 ...
 $ prob                : num  0.469 0.466 0.466 0.466 0.466 ...
 $ score               : num  0.426 0.48 0.481 0.48 0.481 ...
```

- ### Visualization
```
# miRNA-mediated cell-cell communication network
> plot_miRTalk_chord(object = obj)
> plot_miRTalk_circle(object = obj)
> plot_miRTalk_circle_simple(object = obj, celltype = "Tumor")
> plot_miRTalk_heatmap(object = obj)
> plot_miRTalk_sankey(object = obj)
> plot_miR_bubble(object = obj)
> plot_miR_heatmap(object = obj)

# miRNA-target interaction network
> plot_miR2tar_chord(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")
> plot_miR2tar_circle(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")
> plot_miR2tar_heatmap(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")

```
<img src='https://github.com/multitalk/miRTalk/blob/main/img/visualization.png'>

# Note







# About
miRTalk was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn

