# miRTalk
[![R-build-check](https://github.com/shaoxin0801/miRTalk/actions/workflows/R-build-check.yml/badge.svg)](https://github.com/shaoxin0801/miRTalk/actions/workflows/R-build-check.yml) [![version>4](https://img.shields.io/badge/version-%3E%3D4.0-yellow?logo=Rstudio)](#install) [![dockerhub-container](https://img.shields.io/badge/dockerhub-container-yellowgreen?logo=docker)](https://hub.docker.com/repository/docker/zjufanlab/mirtalk)


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
miRTalk method consists of two components, wherein the first is to infer the EV-derived miRNA across cells and the second is to infer the cell-cell communication mediated by EV-derived miRNAs and their downstream targets. Classification and description of miRTalk functions are shown in the __[wiki page](https://github.com/ZJUFanLab/miRTalk/wiki)__
- ### Inference of EV-derived miRNA across cells
```
# object: miRTalk object containg single-cell data matrix and cell types
# mir_info: a data.frame containing the priori knowledge of EV-derived miRNAs

find_miRNA(object, mir_info)
```

- ### Inference of cell-cell communication mediated by EV-derived miRNAs and their downstream targets
```
# object: miRTalk object after running find_miRNA()
# mir2tar: a data.frame containing the priori knowledge of miRNA's target 

find_cci(object, mir2tar)
```

# Note


# About
miRTalk was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn

