# Prerequisites {#prerequisites}

The analysis presented in this book requires a basic understanding of the 
`R` programing language. An introduction to `R` can be found [here](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf) and
in the book [R for Data Science](https://r4ds.hadley.nz/).

Furthermore, it is beneficial to be familiar with single-cell data analysis
using the [Bioconductor](https://www.bioconductor.org/) framework. The 
[Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/) book
gives an excellent overview on data containers and basic analysis that are being
used here.

An overview on IMC as technology and necessary image processing steps can be
found on the [IMC workflow website](https://bodenmillergroup.github.io/IMCWorkflow/). 

Before we get started on IMC data analysis, we will need to make sure that
software dependencies are installed and the example data is downloaded.

## Obtain the code

This book provides R code to perform single-cell and spatial data analysis.
You can copy the individual code chunks into your R scripts or you can obtain
the full code of the book via:

```
git clone https://github.com/BodenmillerGroup/IMCDataAnalysis.git
```

## Software requirements

The R packages needed to execute the presented workflow can either be manually
installed (see section \@ref(manual-install)) or are available within a provided
Docker container (see section \@ref(docker)). The Docker option is useful if you
want to exactly reproduce the presented analysis across operating systems;
however, the manual install gives you more flexibility for exploratory data
analysis.

### Using Docker {#docker}

For reproducibility purposes, we provide a Docker container [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/pkgs/container/imcdataanalysis).

1. After installing [Docker](https://docs.docker.com/get-docker/) you can first pull the container via:

```
docker pull ghcr.io/bodenmillergroup/imcdataanalysis:latest
```

and then run the container:

```
docker run -v /path/to/IMCDataAnalysis:/home/rstudio/IMCDataAnalysis \
	-e PASSWORD=bioc -p 8787:8787  \
	ghcr.io/bodenmillergroup/imcdataanalysis:latest
```

Here, the `/path/to/` needs to be adjusted to where you keep the code and data
of the book.

**Of note: it is recommended to use a date-tagged version of the container to ensure reproducibility**. 
This can be done via:

```
docker pull ghcr.io/bodenmillergroup/imcdataanalysis:<year-month-date>
```

2. An RStudio server session can be accessed via a browser at `localhost:8787` using `Username: rstudio` and `Password: bioc`.  
3. Navigate to `IMCDataAnalysis` and open the `IMCDataAnalysis.Rproj` file.  
4. Code in the individual files can now be executed or the whole workflow can be build by entering `bookdown::render_book()`.

### Manual installation {#manual-install}

The following section describes how to manually install all needed R packages
when not using the provided Docker container.
To install all R packages needed for the analysis, please run:


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rmarkdown", "bookdown", "pheatmap", "viridis", "zoo", 
                       "devtools", "testthat", "tiff", "distill", "ggrepel", 
                       "patchwork", "mclust", "RColorBrewer", "uwot", "Rtsne", 
                       "harmony", "Seurat", "SeuratObject", "cowplot", "kohonen", 
                       "caret", "randomForest", "ggridges", "cowplot", 
                       "gridGraphics", "scales", "tiff", "harmony", "Matrix", 
                       "CATALYST", "scuttle", "scater", "dittoSeq", 
                       "tidyverse", "BiocStyle", "batchelor", "bluster", "scran", 
                       "lisaClust", "spicyR", "iSEE", "imcRtools", "cytomapper",
                       "imcdatasets", "cytoviewer"))

# Github dependencies
devtools::install_github("i-cyto/Rphenograph")
```



### Major package versions

Throughout the analysis, we rely on different R software packages.
This section lists the most commonly used packages in this workflow.

Data containers:

* [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) version 1.12.0
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) version 1.24.0

Data analysis:

* [CATALYST](https://bioconductor.org/packages/release/bioc/html/CATALYST.html) version 1.26.0
* [imcRtools](https://bioconductor.org/packages/release/bioc/html/imcRtools.html) version 1.8.0
* [scuttle](https://bioconductor.org/packages/release/bioc/html/scuttle.html) version 1.12.0
* [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) version 1.30.1
* [batchelor](https://www.bioconductor.org/packages/release/bioc/html/batchelor.html) version 1.18.1
* [bluster](https://www.bioconductor.org/packages/release/bioc/html/bluster.html) version 1.12.0
* [scran](https://www.bioconductor.org/packages/release/bioc/html/scran.html) version 1.30.0
* [harmony](https://github.com/immunogenomics/harmony) version 1.2.0
* [Seurat](https://satijalab.org/seurat/index.html) version 5.0.1
* [lisaClust](https://www.bioconductor.org/packages/release/bioc/html/lisaClust.html) version 1.10.1
* [caret](https://topepo.github.io/caret/) version 6.0.94

Data visualization:

* [cytomapper](https://bioconductor.org/packages/release/bioc/html/cytomapper.html) version 1.14.0
* [cytoviewer](https://bioconductor.org/packages/release/bioc/html/cytoviewer.html) version 1.2.0
* [dittoSeq](https://bioconductor.org/packages/release/bioc/html/dittoSeq.html) version 1.14.0

Tidy R:

* [tidyverse](https://www.tidyverse.org/) version 2.0.0

## Image processing {#image-processing}

The analysis presented here fully relies on packages written in the programming
language `R` and primarily focuses on analysis approaches downstream of image
processing. The example data available at
[https://zenodo.org/record/7575859](https://zenodo.org/record/7575859) were
processed (file type conversion, image segmentation, feature extraction as
explained in Section \@ref(processing)) using the
[steinbock](https://bodenmillergroup.github.io/steinbock/latest/) toolkit. The
exact command line interface calls to process the raw data are shown below:




```bash
#!/usr/bin/env bash
BASEDIR=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "${BASEDIR}"

# raw data collection
mkdir raw
wget https://zenodo.org/record/6449127/files/IMCWorkflow.ilp
wget https://zenodo.org/record/6449127/files/analysis.zip
unzip analysis.zip
rm analysis.zip
rm -r analysis/cpinp
rm -r analysis/cpout
rm -r analysis/histocat
rm -r analysis/ilastik
rm -r analysis/ometiff
cd raw
wget https://zenodo.org/record/5949116/files/panel.csv
wget https://zenodo.org/record/5949116/files/Patient1.zip
wget https://zenodo.org/record/5949116/files/Patient2.zip
wget https://zenodo.org/record/5949116/files/Patient3.zip
wget https://zenodo.org/record/5949116/files/Patient4.zip
cd ${BASEDIR}

# steinbock alias setup
shopt -s expand_aliases
alias steinbock="docker run -v ${BASEDIR}:/data -u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.16.0"

# raw data preprocessing
steinbock preprocess imc panel --namecol Clean_Target
steinbock preprocess imc images --hpf 50

# random forest-based segmentation using Ilastik/CellProfiler
steinbock classify ilastik prepare --cropsize 500 --seed 123
rm pixel_classifier.ilp && mv IMCWorkflow.ilp pixel_classifier.ilp
rm -r ilastik_crops && mv analysis/crops ilastik_crops
steinbock classify ilastik fix --no-backup
steinbock classify ilastik run
steinbock segment cellprofiler prepare
steinbock segment cellprofiler run -o masks_ilastik

# deep learning-based whole-cell segmentation using DeepCell/Mesmer
steinbock segment deepcell --app mesmer --minmax -o masks_deepcell

# single-cell feature extraction
steinbock measure intensities --masks masks_deepcell
steinbock measure regionprops --masks masks_deepcell
steinbock measure neighbors --masks masks_deepcell --type expansion --dmax 4

# data export
steinbock export ome
steinbock export histocat --masks masks_deepcell
steinbock export csv intensities regionprops -o cells.csv
steinbock export csv intensities regionprops --no-concat -o cells_csv
steinbock export fcs intensities regionprops -o cells.fcs
steinbock export fcs intensities regionprops --no-concat -o cells_fcs
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors --no-concat -o cells_h5ad
steinbock export graphs --data intensities

# archiving
zip -r img.zip img
zip -r ilastik_img.zip ilastik_img
zip -r ilastik_crops.zip ilastik_crops
zip -r ilastik_probabilities.zip ilastik_probabilities
zip -r masks_ilastik.zip masks_ilastik
zip -r masks_deepcell.zip masks_deepcell
zip -r intensities.zip intensities
zip -r regionprops.zip regionprops
zip -r neighbors.zip neighbors
zip -r ome.zip ome
zip -r histocat.zip histocat
zip -r cells_csv.zip cells_csv
zip -r cells_fcs.zip cells_fcs
zip -r cells_h5ad.zip cells_h5ad
zip -r graphs.zip graphs
```

## Download example data {#download-data}

Throughout this tutorial, we will access a number of different data types. 
To declutter the analysis scripts, we will already download all needed data here.

To highlight the basic steps of IMC data analysis, we provide example data that
were acquired as part of the **I**ntegrated i**MMU**noprofiling of large adaptive
**CAN**cer patient cohorts projects ([immucan.eu](https://immucan.eu/)). The
raw data of 4 patients can be accessed online at 
[zenodo.org/record/7575859](https://zenodo.org/record/7575859). We will only
download the sample/patient metadata information here:


```r
download.file("https://zenodo.org/record/7575859/files/sample_metadata.csv", 
         destfile = "data/sample_metadata.csv")
```

### Processed multiplexed imaging data

The IMC raw data was either processed using the 
[steinbock](https://github.com/BodenmillerGroup/steinbock) toolkit or the
[IMC Segmentation Pipeline](https://github.com/BodenmillerGroup/ImcSegmentationPipeline).
Image processing included file type conversion, cell segmentation and feature
extraction. 

**steinbock output**

This book uses the output of the `steinbock` framework when applied to process
the example data. The processed data includes the single-cell mean intensity
files, the single-cell morphological features and spatial locations, spatial
object graphs in form of edge lists indicating cells in close proximity, hot
pixel filtered multi-channel images, segmentation masks, image metadata and
channel metadata. All these files will be downloaded here for later use. The
commands which were used to generate this data can be found in the shell script
above.


```r
# download intensities
url <- "https://zenodo.org/record/7624451/files/intensities.zip"
destfile <- "data/steinbock/intensities.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/steinbock", overwrite=TRUE)
unlink(destfile)

# download regionprops
url <- "https://zenodo.org/record/7624451/files/regionprops.zip"
destfile <- "data/steinbock/regionprops.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/steinbock", overwrite=TRUE)
unlink(destfile)

# download neighbors
url <- "https://zenodo.org/record/7624451/files/neighbors.zip"
destfile <- "data/steinbock/neighbors.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/steinbock", overwrite=TRUE)
unlink(destfile)

# download images
url <- "https://zenodo.org/record/7624451/files/img.zip"
destfile <- "data/steinbock/img.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/steinbock", overwrite=TRUE)
unlink(destfile)

# download masks
url <- "https://zenodo.org/record/7624451/files/masks_deepcell.zip"
destfile <- "data/steinbock/masks_deepcell.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/steinbock", overwrite=TRUE)
unlink(destfile)

# download individual files
download.file("https://zenodo.org/record/7624451/files/panel.csv", 
              "data/steinbock/panel.csv")
download.file("https://zenodo.org/record/7624451/files/images.csv", 
              "data/steinbock/images.csv")
download.file("https://zenodo.org/record/7624451/files/steinbock.sh", 
              "data/steinbock/steinbock.sh")
```

**IMC Segmentation Pipeline output**

The example data was also processed using the 
[IMC Segmetation Pipeline](https://github.com/BodenmillerGroup/ImcSegmentationPipeline) (version 3). 
To highlight the use of the reader function for this type of output, we will need
to download the `cpout` folder which is part of the `analysis` folder. The `cpout`
folder stores all relevant output files of the pipeline. For a full description
of the pipeline, please refer to the [docs](https://bodenmillergroup.github.io/ImcSegmentationPipeline/).


```r
# download analysis folder
url <- "https://zenodo.org/record/7997296/files/analysis.zip"
destfile <- "data/ImcSegmentationPipeline/analysis.zip"
download.file(url, destfile)
unzip(destfile, exdir="data/ImcSegmentationPipeline", overwrite=TRUE)
unlink(destfile)

unlink("data/ImcSegmentationPipeline/analysis/cpinp/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/crops/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/histocat/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/ilastik/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/ometiff/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/cpout/images/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/cpout/probabilities/", recursive=TRUE)
unlink("data/ImcSegmentationPipeline/analysis/cpout/masks/", recursive=TRUE)
```

### Files for spillover matrix estimation

To highlight the estimation and correction of channel-spillover as described by
[@Chevrier2017], we can access an example spillover-acquisition from:


```r
download.file("https://zenodo.org/record/7575859/files/compensation.zip",
              "data/compensation.zip")
unzip("data/compensation.zip", exdir="data", overwrite=TRUE)
unlink("data/compensation.zip")
```

### Gated cells

In Section \@ref(classification), we present a cell type classification approach
that relies on previously gated cells. This ground truth data is available
online at [zenodo.org/record/8095133](https://zenodo.org/record/8095133) and
will be downloaded here for later use:


```r
download.file("https://zenodo.org/record/8095133/files/gated_cells.zip",
              "data/gated_cells.zip")
unzip("data/gated_cells.zip", exdir="data", overwrite=TRUE)
unlink("data/gated_cells.zip")
```

## Software versions {#sessionInfo}

<details>
   <summary>SessionInfo</summary>
   

```
## R version 4.3.2 (2023-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 22.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Etc/UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] cytoviewer_1.2.0            caret_6.0-94               
##  [3] lattice_0.21-9              lisaClust_1.10.1           
##  [5] scran_1.30.0                bluster_1.12.0             
##  [7] lubridate_1.9.3             forcats_1.0.0              
##  [9] stringr_1.5.1               dplyr_1.1.4                
## [11] purrr_1.0.2                 readr_2.1.4                
## [13] tidyr_1.3.0                 tibble_3.2.1               
## [15] tidyverse_2.0.0             dittoSeq_1.14.0            
## [17] cytomapper_1.14.0           EBImage_4.44.0             
## [19] imcRtools_1.8.0             scater_1.30.1              
## [21] ggplot2_3.4.4               scuttle_1.12.0             
## [23] SpatialExperiment_1.12.0    CATALYST_1.26.0            
## [25] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
## [27] Biobase_2.62.0              GenomicRanges_1.54.1       
## [29] GenomeInfoDb_1.38.5         IRanges_2.36.0             
## [31] S4Vectors_0.40.2            BiocGenerics_0.48.1        
## [33] MatrixGenerics_1.14.0       matrixStats_1.2.0          
## 
## loaded via a namespace (and not attached):
##   [1] vroom_1.6.5                 tiff_0.1-12                
##   [3] nnet_7.3-19                 goftest_1.2-3              
##   [5] DT_0.31                     HDF5Array_1.30.0           
##   [7] TH.data_1.1-2               vctrs_0.6.5                
##   [9] spatstat.random_3.2-2       digest_0.6.33              
##  [11] png_0.1-8                   shape_1.4.6                
##  [13] proxy_0.4-27                ggrepel_0.9.4              
##  [15] spicyR_1.14.2               deldir_2.0-2               
##  [17] parallelly_1.36.0           magick_2.8.2               
##  [19] MASS_7.3-60                 reshape2_1.4.4             
##  [21] httpuv_1.6.13               foreach_1.5.2              
##  [23] withr_2.5.2                 xfun_0.41                  
##  [25] ggpubr_0.6.0                ellipsis_0.3.2             
##  [27] survival_3.5-7              RTriangle_1.6-0.12         
##  [29] ggbeeswarm_0.7.2            RProtoBufLib_2.14.0        
##  [31] drc_3.0-1                   systemfonts_1.0.5          
##  [33] zoo_1.8-12                  GlobalOptions_0.1.2        
##  [35] gtools_3.9.5                promises_1.2.1             
##  [37] rstatix_0.7.2               globals_0.16.2             
##  [39] rhdf5filters_1.14.1         rhdf5_2.46.1               
##  [41] miniUI_0.1.1.1              archive_1.1.7              
##  [43] units_0.8-5                 generics_0.1.3             
##  [45] concaveman_1.1.0            zlibbioc_1.48.0            
##  [47] ScaledMatrix_1.10.0         ggraph_2.1.0               
##  [49] polyclip_1.10-6             GenomeInfoDbData_1.2.11    
##  [51] SparseArray_1.2.3           fftwtools_0.9-11           
##  [53] xtable_1.8-4                doParallel_1.0.17          
##  [55] evaluate_0.23               S4Arrays_1.2.0             
##  [57] hms_1.1.3                   bookdown_0.37              
##  [59] irlba_2.3.5.1               colorspace_2.1-0           
##  [61] spatstat.data_3.0-3         magrittr_2.0.3             
##  [63] later_1.3.2                 viridis_0.6.4              
##  [65] spatstat.geom_3.2-7         future.apply_1.11.1        
##  [67] XML_3.99-0.16               cowplot_1.1.2              
##  [69] class_7.3-22                svgPanZoom_0.3.4           
##  [71] pillar_1.9.0                nlme_3.1-163               
##  [73] iterators_1.0.14            compiler_4.3.2             
##  [75] beachmat_2.18.0             shinycssloaders_1.0.0      
##  [77] stringi_1.8.3               gower_1.0.1                
##  [79] sf_1.0-15                   tensor_1.5                 
##  [81] minqa_1.2.6                 ClassifyR_3.6.2            
##  [83] plyr_1.8.9                  crayon_1.5.2               
##  [85] abind_1.4-5                 locfit_1.5-9.8             
##  [87] sp_2.1-2                    graphlayouts_1.0.2         
##  [89] bit_4.0.5                   terra_1.7-65               
##  [91] sandwich_3.1-0              codetools_0.2-19           
##  [93] multcomp_1.4-25             recipes_1.0.9              
##  [95] BiocSingular_1.18.0         bslib_0.6.1                
##  [97] e1071_1.7-14                GetoptLong_1.0.5           
##  [99] mime_0.12                   MultiAssayExperiment_1.28.0
## [101] splines_4.3.2               circlize_0.4.15            
## [103] Rcpp_1.0.11                 sparseMatrixStats_1.14.0   
## [105] knitr_1.45                  utf8_1.2.4                 
## [107] clue_0.3-65                 lme4_1.1-35.1              
## [109] listenv_0.9.0               nnls_1.5                   
## [111] DelayedMatrixStats_1.24.0   ggsignif_0.6.4             
## [113] Matrix_1.6-4                scam_1.2-14                
## [115] statmod_1.5.0               tzdb_0.4.0                 
## [117] svglite_2.1.3               tweenr_2.0.2               
## [119] pkgconfig_2.0.3             pheatmap_1.0.12            
## [121] tools_4.3.2                 cachem_1.0.8               
## [123] viridisLite_0.4.2           DBI_1.2.0                  
## [125] numDeriv_2016.8-1.1         fastmap_1.1.1              
## [127] rmarkdown_2.25              scales_1.3.0               
## [129] grid_4.3.2                  shinydashboard_0.7.2       
## [131] broom_1.0.5                 sass_0.4.8                 
## [133] carData_3.0-5               rpart_4.1.21               
## [135] farver_2.1.1                tidygraph_1.3.0            
## [137] mgcv_1.9-0                  yaml_2.3.8                 
## [139] cli_3.6.2                   lifecycle_1.0.4            
## [141] mvtnorm_1.2-4               lava_1.7.3                 
## [143] backports_1.4.1             BiocParallel_1.36.0        
## [145] cytolib_2.14.0              timechange_0.2.0           
## [147] gtable_0.3.4                rjson_0.2.21               
## [149] ggridges_0.5.5              parallel_4.3.2             
## [151] pROC_1.18.5                 limma_3.58.1               
## [153] colourpicker_1.3.0          jsonlite_1.8.8             
## [155] edgeR_4.0.3                 bitops_1.0-7               
## [157] bit64_4.0.5                 Rtsne_0.17                 
## [159] FlowSOM_2.10.0              spatstat.utils_3.0-4       
## [161] BiocNeighbors_1.20.1        flowCore_2.14.0            
## [163] jquerylib_0.1.4             metapod_1.10.1             
## [165] dqrng_0.3.2                 timeDate_4032.109          
## [167] shiny_1.8.0                 ConsensusClusterPlus_1.66.0
## [169] htmltools_0.5.7             distances_0.1.10           
## [171] glue_1.6.2                  XVector_0.42.0             
## [173] RCurl_1.98-1.13             classInt_0.4-10            
## [175] jpeg_0.1-10                 gridExtra_2.3              
## [177] boot_1.3-28.1               igraph_1.6.0               
## [179] R6_2.5.1                    cluster_2.1.4              
## [181] Rhdf5lib_1.24.1             ipred_0.9-14               
## [183] nloptr_2.0.3                DelayedArray_0.28.0        
## [185] tidyselect_1.2.0            vipor_0.4.7                
## [187] plotrix_3.8-4               ggforce_0.4.1              
## [189] raster_3.6-26               car_3.1-2                  
## [191] future_1.33.1               ModelMetrics_1.2.2.2       
## [193] rsvd_1.0.5                  munsell_0.5.0              
## [195] KernSmooth_2.23-22          data.table_1.14.10         
## [197] htmlwidgets_1.6.4           ComplexHeatmap_2.18.0      
## [199] RColorBrewer_1.1-3          rlang_1.1.2                
## [201] spatstat.sparse_3.0-3       spatstat.explore_3.2-5     
## [203] lmerTest_3.1-3              colorRamps_2.3.1           
## [205] ggnewscale_0.4.9            fansi_1.0.6                
## [207] hardhat_1.3.0               beeswarm_0.4.0             
## [209] prodlim_2023.08.28
```
</details>



