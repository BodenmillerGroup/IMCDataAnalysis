# Docker inheritance
FROM rocker/rstudio:latest

RUN apt-get -y update \
    && apt-get install -y --no-install-recommends apt-utils \
    && apt-get install -y --no-install-recommends zlib1g-dev libglpk-dev libmagick++-dev libfftw3-dev libxml2-dev libxt-dev curl libcairo2-dev libproj-dev libgdal-dev libudunits2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/ilists/*

RUN R -e 'install.packages(c("rmarkdown", "bookdown", "pheatmap", "viridis", "zoo", "BiocManager", "devtools", "testthat", "tiff", \
                             "distill", "ggrepel", "patchwork", "mclust", "RColorBrewer", "uwot", "Rtsne", "harmony", \
                             "Seurat", "SeuratObject", "cowplot", "kohonen", "caret", "randomForest", "ggridges", "cowplot", \
                             "gridGraphics", "scales", "tiff", "harmony", "Matrix"))'
RUN R -e 'BiocManager::install(c("CATALYST", "scuttle", "scater", "dittoSeq", "tidyverse", "BiocStyle", "batchelor", "bluster", \
                                 "scran", "lisaClust", "spicyR", "iSEE", "imcRtools", "cytomapper", "imcdatasets", "cytoviewer"")'
RUN R -e 'devtools::install_github(c("i-cyto/Rphenograph"))'

