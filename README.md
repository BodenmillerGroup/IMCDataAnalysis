## An end-to-end workflow for multiplexed image processing and analysis

This branch of the repository contains the code to reproduce the analysis presented in the following paper:

```
TODO upon acceptance
```

The following lists the system and software requirements to run the workflow.

### System requirements

To run the workflow, a computer with a recent version of a Windows, Mac, or Linux operating system (OS) is required. 
With increasing dataset size, more memory is required and we recommend at least 8 GB RAM. 
Alternatively, a high performance computer (e.g. cluster) can be used, provided Docker can be installed (see below). 
For this manuscript, the workflow was run on MacOS Catalina (10.15.7), 2.7 GHz Quad-Core Intel Core i7, 16 GB 2133 MHz LPDDR3.

### Software used

* **napari & napari-imc (IMC-specific):** The multi-dimensional image viewer napari  (https://napari.org) together with the napari-imc plugin for loading imaging mass cytometry files (https://github.com/BodenmillerGroup/napari-imc) were used to visualize and inspect raw multiplexed imaging data. Python 3.9.12 (https://www.python.org), napari 0.4.16, and napari-imc 0.6.5 were installed into a fresh conda (https://conda.io) environment; see below for installation instructions.
* **steinbock Docker container:** The multi-channel image processing toolkit steinbock  (https://bodenmillergroup.github.io/steinbock) was used to pre-process multiplexed imaging data, perform image segmentation, and extract single-cell data. The steinbock Docker container v0.14.1 was pulled from the GitHub container registry using Docker Desktop 4.9.0 for Mac; see below for installation instructions.
* **Ilastik/CellProfiler-based segmentation pipeline:** Multiplexed image processing using random forest-based pixel classification and watershed-based cell segmentation was performed using the Ilastik/CellProfiler-based segmentation pipeline v3.4 (https://bodenmillergroup.github.io/ImcSegmentationPipeline/); see below for installation instructions.

In addition, in order to use the pipeline, the following software need to be installed:
* **Ilastik:** The Ilastik software24 is used for pixel-classification prior to cell segmentation and can be installed from https://www.ilastik.org/download.html. The version used for this workflow is v1.3.3post3.
* **CellProfiler:** The CellProfiler software25 is used to segment individual cells. The tool can be installed from https://cellprofiler.org/previous-releases on Windows (64-bit) and MacOS (10.14+). The version used in this workflow is v4.2.1.
* **R setup:** Downstream analysis after image processing is conducted using the statistical programming language R, which can be installed from https://cran.r-project.org/ following the OS-specific instructions. The version used in this workflow is v4.2.1. 
* The RStudio software offers an easy-to-use GUI for data analysis in R. It can be installed from https://www.rstudio.com/products/rstudio/download/.
* The versions of all R libraries used in this workflow can be seen below. 

```
TODO add session info
```

### Installation instructions

* **napari & napari-imc:** Install the conda package manager according to the instructions at https://docs.conda.io/projects/conda/en/latest/user-guide/install/ 
Create a new conda environment with Python 3.9:
```
conda create -n napari-imc python=3.9
```
Activate the conda environment and Install napari & napari-imc:
```
conda activate napari-imc
pip install “napari[all]==0.4.16” napari-imc==0.6.5
```
* **steinbock:** Instructions to install the dockerized steinbock toolkit can be found at https://bodenmillergroup.github.io/steinbock/v0.14.1/install-docker/. In particular, to run the steinbock container, Docker needs to be installed first (see online instructions). For this manuscript, we run steinbock using the following alias:
```
alias steinbock="docker run -v /path/to/data/steinbock:/data -u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.14.1"
```
CRITICAL: In the command above the `/path/to/data/steinbock` needs to be adapted and replaced by the anticipated working directory. 

* **Ilastik/CellProfiler-based segmentation pipeline:** the pre-processing steps of the pipeline are performed in Python using a custom script. To setup the pre-processing script, the following steps need to be performed:

Install conda from https://docs.conda.io/projects/conda/en/latest/user-guide/install/

Clone the repository
```
git clone --recursive      https://github.com/BodenmillerGroup/ImcSegmentationPipeline.git
```
Setup the imcsegpipe conda environment:
```
cd ImcSegmentationPipeline
conda env create -f environment.yml
```

Configure CellProfiler to use the required plugins by opening the CellProfiler GUI, selecting Preferences and setting the CellProfiler plugins directory to path/to/ImcSegmentationPipeline/resources/ImcPluginsCP/plugins and restart CellProfiler. 

* **R libraries:** the workflow highlights the use of the cytomapper and imcRtools R/Bioconductor packages. With the following callcall below, we can install the development versions of both packages can be installed:.

```
if (!requireNamespace("devtools", quietly = TRUE))
install.packages("devtools")

devtools::install_github(c("BodenmillerGroup/imcRtools", 
                           "BodenmillerGroup/cytomapper"))
```

In addition, a number of additional R packages need to be installed to follow the workflow:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install(c("pheatmap", "viridis",
                       "tiff", "distill", "openxlsx", "ggrepel", "patchwork",
                       "mclust", "RColorBrewer", "uwot", "Rtsne", "caret",                                                
                       "randomForest", "ggridges", "gridGraphics", "scales", 
                       "CATALYST", "scuttle", "scater", "dittoSeq", 
                       "tidyverse", "batchelor", "bluster","scran"))
```

To install the required software around 1-2 hours need to be taken into account.

### Running the workflow

For reproducing the exact analysis, we recommend to have R v4.2.1 installed and Biocondcutor release version 3.15. 
In RStudio open the `protocol.Rmd` file, start the docker engine and proceed with code execution.
Applying the workflow to the provided dataset takes roughly 30 minutes and provides the raw data files, data generated by the `steinbock` toolkit and a `SpatialExperiment` object storing all analysis results.
