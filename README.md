# miRTalk
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen?logo=github)](https://github.com/multitalk/miRTalk/actions)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13856217.svg)](https://doi.org/10.5281/zenodo.13856217)

### Extracellular vesicle-derived miRNA-mediated cell-cell communication inference for single-cell transcriptomic data

<img src='https://github.com/multitalk/miRTalk/blob/main/img/workflow.png'>

MicroRNAs are released from cells in extracellular vesicles (EVs), including exosomes and microvesicles, representing an essential mode of cell-cell communications via regulatory effect on gene expression of receivers. The advent of single-cell RNA-sequencing (scRNA-seq) technologies has ushered in an era of elucidating EV-derived miRNA-mediated cell-cell communications at an unprecedented resolution. However, the lack of computational method to infer such communications from scRNA-seq data poses an outstanding challenge. Herein, we present miRTalk, a pioneering framework for inferreing EV-derived miRNA-mediated cell-cell communications with a curated database, miRTalkDB, which include EV-derived miRNA-target associations. Our miRTalk considers 1) the potential of producing EVs and 2) the expression of miRNAs in senders as well as 3) the activation of miRNA processing machinery and 4) the expression of target genes in receivers.


# Install

```
# install devtools and install
install.packages(pkgs = "devtools"")
devtools::install_github("multitalk/miRTalk")
```

OR

```
# download the repository as ZIP
devtools::install_local("/path/to/miRTalk-main.zip")
```

# Usage
miRTalk method consists of three components, wherein the first is to infer the EV-derived miRNA across cells and the highly variable target genes, the second is to infer the cell-cell communication mediated by EV-derived miRNAs and their downstream targets. The third part is to visualize the miRNA-mediated cell-cell communication network and miRNA-target interaction network. Please refer to the __[tutorial vignette](https://raw.githack.com/multitalk/miRTalk/main/vignettes/tutorial.html)__ with demo data processing steps. Detailed functions see the __[document](https://raw.githack.com/multitalk/miRTalk/main/vignettes/miRTalk.pdf)__.Step-by-step procedures can be found under the __[wiki page](https://github.com/multitalk/miRTalk/wiki)__.

<img src='https://github.com/multitalk/miRTalk/blob/main/img/visualization.png'>

- ### Inference of EV-derived miRNA and highly variable target genes
```
# sc_data: a data.frame or matrix or dgCMatrix containing raw counts of single-cell RNA-seq data
# sc_celltype: a character containing the cell type of the single-cell RNA-seq data
# condition: a character with the same length as the number of cells, e.g., control/disease/treatment, phase 1/2/3, men/women
# evbiog: a data.frame of the system data containing extracellular vesicle biogenesis genes (GO:0140112)
# risc: a data.frame of the system data containing RNA-induced silencing complex (RISC) related genes
# mir_info: 
# 
> obj <- create_miRTalk(sc_data = sc_data,
                         sc_celltype = sc_celltype,
                         species = "Human",
                         condition = rep("cancer", length(sc_celltype)),
                         evbiog = evbiog,risc = risc)
> obj
An object of class miRTalk
0 EV-derived miRNA-target interactions

> obj <- find_hvtg(object = obj)
> obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
```

- ### Inference of cell-cell communication mediated by EV-derived miRNA-target interactions
```
# object: miRTalk object after running find_hvtg() abd find_miRNA() 

> obj <- find_miRTalk(obj, if_doParallel = F)
[cancer] 
[++++++++++++++++++++++++++++++] Finished:100% time:00:11:00
> obj
An object of class miRTalk 
460 EV-derived miRNA-target interactions

> obj_cci <- get_miRTalk_cci(obj)
> str(obj_cci)
'data.frame':    449 obs. of  10 variables:
 $ celltype_sender     : chr  "Bcell" "Myeloid" "Tcell" "Tumor" ...
 $ celltype_receiver   : chr  "Bcell" "Bcell" "Bcell" "Bcell" ...
 $ miRNA               : chr  "hsa-miR-3916" "hsa-miR-3916" "hsa-miR-3916" "hsa-miR-3916" ...
 $ EVmiR_score         : num  0.0174 0.0238 0.0313 0.0448 0.0117 ...
 $ target_gene         : chr  "CLDN4" "CLDN4" "CLDN4" "CLDN4" ...
 $ target_gene_activity: num  0.00086 0.00086 0.00086 0.00086 0.00989 ...
 $ score               : num  0.01266 0.04056 0.03305 0.01867 0.00402 ...
 $ condition           : chr  "condition" "condition" "condition" "condition" ...
 $ miR2tar             : chr  "hsa-miR-3916:CLDN4" "hsa-miR-3916:CLDN4" "hsa-miR-3916:CLDN4" "hsa-miR-3916:CLDN4" ...
 $ specifity           : num  0.148 0.203 0.267 0.382 0.246 ...

- ### Visualization of miRNA-mediated cell-cell communication network

```
> plot_miRTalk_chord(object = obj)
> plot_miRTalk_circle(object = obj)
> plot_miRTalk_circle_simple(object = obj, celltype = "Tumor")
> plot_miRTalk_heatmap(object = obj)
> plot_miRTalk_sankey(object = obj)
> plot_miR_bubble(object = obj)
> plot_miR_heatmap(object = obj)
```

- ### Visualization of specifity analysis

```
> plot_target_heatmap(object = obj, celltype = "Bcell")
```

- ### Visualization of miRNA-target interaction network
```
> plot_miR2tar_chord(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")
> plot_miR2tar_circle(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")
> plot_miR2tar_heatmap(object = obj, celltype_sender = "Tumor", celltype_receiver = "Stromal")
```

# Note
[![miRTalkDB](https://img.shields.io/badge/miRTalkDB-v1.0-yellow)](https://github.com/multitalk/miRTalk/tree/main/data) [![miRTalk-tutorial](https://img.shields.io/badge/miRTalk-tutorial-blue)](https://raw.githack.com/multitalk/miRTalk/main/vignettes/tutorial.html)

- __miRTalk can be applied to either [single-cell transcriptomic data](https://github.com/multitalk/miRTalk/tree/main/inst/extdata) or [spatial transcriptomic data](https://doi.org/10.1038/s41587-022-01517-6)__
- __miRTalk allows to infer miRNA-target interactions with [positive regulation](https://github.com/multitalk/miRTalk/wiki/Inference-of-miRNA%E2%80%90target-interactions-with-positive-regulation). The default is [negative regulation](https://github.com/multitalk/miRTalk/wiki/Inference-of-miRNA%E2%80%90target-interactions-with-negative-regulation)__
- miRTalk allows to use the parallel processing for `find_miRTalk`
- __miRTalk allows to infer miRNA-target interactions [under different conditions](https://github.com/multitalk/miRTalk/wiki/Inference-of-miRNA%E2%80%90target-interactions-under-different-conditions)__
- __miRTalk allows to use [different databases of miRNA-target interactions](https://github.com/multitalk/miRTalk/wiki/Inference-with-different-databases-of-miRNA%E2%80%90target-interactions)__
- __miRTalk allows to infer miRNA-target interactions [with custom data provided by users](https://github.com/multitalk/miRTalk/wiki/Inference-with-custom-data-provided-by-users)__
  - custom evbiog, risc, mir_data, mir2tar, highly varible genes, gene2path, mir2path
- __miRTalk allows to infer miRNA-target interactions [with bulk RNA-seq data](https://github.com/multitalk/miRTalk/wiki/Inference-of-miRNA%E2%80%90target-interactions-with-bulk-RNA%E2%80%90seq-data)__
- __miRTalk allows to [use EXOmotif](https://github.com/multitalk/miRTalk/wiki/Use-EXOmotif-to-refine-miRNA%E2%80%90target-interactions) to refine miRNA-target interactions__
- __miRTalk allows to [use human data for the inference of mouse and rat data](https://github.com/multitalk/miRTalk/wiki/Use-human-data-for-the-inference-of-mouse-and-rat-data)__
- __miRTalk can provide [analysis of specifity for a given receiver cell type](https://github.com/multitalk/miRTalk/wiki/Analysis-of-specifity-for-a-given-receiver-cell-type)__
- __miRTalk can provide [analysis of functional annotations for a miRNA and its target gene](https://github.com/multitalk/miRTalk/wiki/Analysis-of-functional-annotations-for-a-miRNA-and-its-target-gene)__
- __miRTalk can provide [analysis of potential of circulating miRNAs and organ-organ communication](https://github.com/multitalk/miRTalk/wiki/Analysis-of-potential-of-circulating-miRNAs-and-organ%E2%80%90organ-communication)__
- __miRTalkDB can be download at[`data/`](https://github.com/multitalk/miRTalk/tree/main/data)__

__Please refer to the [tutorial vignette](https://raw.githack.com/multitalk/miRTalk/main/vignettes/tutorial.html) with demo data processing steps. Step-by-step procedures can be found under the [wiki page](https://github.com/multitalk/miRTalk/wiki). Detailed functions see the [document](https://raw.githack.com/multitalk/miRTalk/main/vignettes/miRTalk.pdf)__

# About
miRTalk was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn

