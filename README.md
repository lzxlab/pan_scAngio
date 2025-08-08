# pan_scAngio
The central code for the paper "Pancancer pro-angiogenic atlas unravels tumor-educated pericyte-augmented anti-angiogenic resistance"

# Overview

Targeting angiogenesis represents a cornerstone in the development of antitumoral therapy, yet the suboptimal responses and inevitable resistance significantly limit its clinical efficacy. To better understand the mechanisms underlying anti-angiogenic resistance, this study systematically characterized the pancancer angiogenic landscape based on 1.24 million individual cells across 13 cancer types. Our analyses revealed that perivascular cells (PCs) and PC-derived non-classical angiogenic factors PGF/ANGPT2 correlate with angiogenesis more prominently than the classical VEGFA-centric model. Tumor-educated MCAM+ immature PCs (imPCs) were spatio-transcriptionally identified as the primary source of PGF/ANGPT2, driving alternative angiogenesis and serving as a major contributor to ɑVEGFR resistance. To address the limitations associated with single-targeted endothelial cells (ECs) inhibition via ɑVEGFR, we developed an innovative dual-targeting strategy that combines ɑVEGFR with MCAM-directed antibody-drug conjugates (MCAM-ADCs) to inhibit ECs and PCs simultaneously. This approach demonstrated superior tumor control, offering MCAM-ADC as a promising translational solution to circumvent anti-angiogenic resistance. Overall, our findings redefine the mechanisms of angiogenic resistance and suggest an innovative dual EC/PC inhibition strategy for more effective anticancer therapy.

# Repo Contents

- [R](./R): `R` package code.


# System Requirements

## Hardware Requirements

The code requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 100 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 100+ GB  
CPU: 20+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (100 GB RAM, 20 cores@3.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Centos 8  
Mac OSX:  
Windows:  

The CRAN package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the package, users should have `R` version 4.2.0 or higher, and several packages set up from CRAN.

#### Installing R version 4.2.0 on Centos 8 

the latest version of R can be installed by using Conda:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
conda --version
conda create --name r_42 r-base=4.2.0 r-essentials
```

which should install in about 20 seconds.

# Installation Guide


### Package dependencies

Users should install the following packages prior to installing `lolR`, from an `R` terminal:

```
conda activate r_42
R
if (!require("BiocManager", quietly = TRUE)){
install.packages("BiocManager")
}

BiocManager::install(c("Seurat","data.table","ggplot2","ggridges","dplyr","cowplot","PCAtools","ggpubr","devtools","monocle","monocle3"))
BiocManager::install(c("Japrin/sscVis","Japrin/sscClust","Japrin/scPip"))
```
which will install in about 30 minutes on a machine with the recommended specs.

### Package Installation and Using the packge:
```
git clone lzxlab/pan_scAngio
Rscript Figure_1.R
```
Users can check [Github snapshot](https://github.com/lzxlab/pan_scAngio/) for details. The versions of software are, specifically:

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/lzxlab/pan_scAngio/issues). 
