## An end-to-end workflow for multiplexed image processing and analysis

This folder of the repository contains the code to reproduce the analysis presented in the following paper:

```
TODO upon acceptance
```

The following lists the system and software requirements to run the workflow.

### Reproducing the analysis



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
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] scales_1.2.1                pheatmap_1.0.12            
##  [3] igraph_1.3.5                viridis_0.6.2              
##  [5] viridisLite_0.4.1           caret_6.0-93               
##  [7] lattice_0.20-45             scran_1.24.1               
##  [9] BiocParallel_1.30.3         bluster_1.6.0              
## [11] batchelor_1.12.3            scater_1.24.0              
## [13] scuttle_1.6.3               patchwork_1.1.2            
## [15] dittoSeq_1.8.1              CATALYST_1.20.1            
## [17] cytomapper_1.9.1            EBImage_4.38.0             
## [19] forcats_0.5.2               stringr_1.4.1              
## [21] dplyr_1.0.10                purrr_0.3.4                
## [23] readr_2.1.2                 tidyr_1.2.1                
## [25] tibble_3.1.8                ggplot2_3.3.6              
## [27] tidyverse_1.3.2             openxlsx_4.2.5             
## [29] imcRtools_1.3.7             SpatialExperiment_1.6.1    
## [31] SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1
## [33] Biobase_2.56.0              GenomicRanges_1.48.0       
## [35] GenomeInfoDb_1.32.4         IRanges_2.30.1             
## [37] S4Vectors_0.34.0            BiocGenerics_0.42.0        
## [39] MatrixGenerics_1.8.1        matrixStats_0.62.0         
## [41] BiocStyle_2.24.0           
## 
## loaded via a namespace (and not attached):
##   [1] rsvd_1.0.5                  svglite_2.1.0              
##   [3] class_7.3-20                fftwtools_0.9-11           
##   [5] V8_4.2.1                    foreach_1.5.2              
##   [7] crayon_1.5.1                MASS_7.3-58.1              
##   [9] rhdf5filters_1.8.0          nlme_3.1-159               
##  [11] backports_1.4.1             reprex_2.0.2               
##  [13] rlang_1.0.6                 XVector_0.36.0             
##  [15] readxl_1.4.1                irlba_2.3.5                
##  [17] limma_3.52.3                rjson_0.2.21               
##  [19] bit64_4.0.5                 glue_1.6.2                 
##  [21] parallel_4.2.1              vipor_0.4.5                
##  [23] classInt_0.4-7              shinydashboard_0.7.2       
##  [25] haven_2.5.1                 tidyselect_1.1.2           
##  [27] distances_0.1.8             XML_3.99-0.10              
##  [29] zoo_1.8-11                  sf_1.0-8                   
##  [31] ggpubr_0.4.0                nnls_1.4                   
##  [33] xtable_1.8-4                magrittr_2.0.3             
##  [35] evaluate_0.16               cli_3.4.1                  
##  [37] zlibbioc_1.42.0             rstudioapi_0.14            
##  [39] sp_1.5-0                    rpart_4.1.16               
##  [41] bslib_0.4.0                 shiny_1.7.2                
##  [43] BiocSingular_1.12.0         xfun_0.33                  
##  [45] clue_0.3-61                 cluster_2.1.4              
##  [47] tidygraph_1.2.2             ggrepel_0.9.1              
##  [49] listenv_0.8.0               future_1.28.0              
##  [51] png_0.1-7                   ipred_0.9-13               
##  [53] withr_2.5.0                 bitops_1.0-7               
##  [55] aws.signature_0.6.0         ggforce_0.3.4              
##  [57] RBGL_1.72.0                 plyr_1.8.7                 
##  [59] cellranger_1.1.0            ncdfFlow_2.42.1            
##  [61] RTriangle_1.6-0.10          hardhat_1.2.0              
##  [63] pROC_1.18.0                 e1071_1.7-11               
##  [65] dqrng_0.3.0                 pillar_1.8.1               
##  [67] RcppParallel_5.1.5          GlobalOptions_0.1.2        
##  [69] cachem_1.0.6                multcomp_1.4-20            
##  [71] fs_1.5.2                    CytoML_2.8.1               
##  [73] raster_3.6-3                GetoptLong_1.0.5           
##  [75] DelayedMatrixStats_1.18.0   vctrs_0.4.1                
##  [77] ellipsis_0.3.2              generics_0.1.3             
##  [79] lava_1.6.10                 tools_4.2.1                
##  [81] beeswarm_0.4.0              munsell_0.5.0              
##  [83] tweenr_2.0.2                aws.s3_0.3.21              
##  [85] proxy_0.4-27                DelayedArray_0.22.0        
##  [87] fastmap_1.1.0               compiler_4.2.1             
##  [89] abind_1.4-5                 httpuv_1.6.6               
##  [91] prodlim_2019.11.13          GenomeInfoDbData_1.2.8     
##  [93] gridExtra_2.3               edgeR_3.38.4               
##  [95] ggnewscale_0.4.7            ggpointdensity_0.1.0       
##  [97] deldir_1.0-6                utf8_1.2.2                 
##  [99] later_1.3.0                 recipes_1.0.1              
## [101] jsonlite_1.8.0              concaveman_1.1.0           
## [103] graph_1.74.0                ScaledMatrix_1.4.1         
## [105] carData_3.0-5               sparseMatrixStats_1.8.0    
## [107] promises_1.2.0.1            car_3.1-0                  
## [109] doParallel_1.0.17           latticeExtra_0.6-30        
## [111] R.utils_2.12.0              rmarkdown_2.16             
## [113] sandwich_3.0-2              cowplot_1.1.1              
## [115] textshaping_0.3.6           statmod_1.4.37             
## [117] Rtsne_0.16                  uwot_0.1.14                
## [119] HDF5Array_1.24.2            survival_3.4-0             
## [121] ResidualMatrix_1.6.1        yaml_2.3.5                 
## [123] plotrix_3.8-2               systemfonts_1.0.4          
## [125] cytolib_2.8.0               flowWorkspace_4.8.0        
## [127] htmltools_0.5.3             locfit_1.5-9.6             
## [129] graphlayouts_0.8.1          digest_0.6.29              
## [131] assertthat_0.2.1            mime_0.12                  
## [133] tiff_0.1-11                 units_0.8-0                
## [135] future.apply_1.9.1          data.table_1.14.2          
## [137] R.oo_1.25.0                 flowCore_2.8.0             
## [139] drc_3.0-1                   ragg_1.2.2                 
## [141] splines_4.2.1               labeling_0.4.2             
## [143] Rhdf5lib_1.18.2             googledrive_2.0.0          
## [145] RCurl_1.98-1.8              broom_1.0.1                
## [147] hms_1.1.2                   modelr_0.1.9               
## [149] rhdf5_2.40.0                colorspace_2.0-3           
## [151] DropletUtils_1.16.0         ConsensusClusterPlus_1.60.0
## [153] base64enc_0.1-3             BiocManager_1.30.18        
## [155] ggbeeswarm_0.6.0            shape_1.4.6                
## [157] nnet_7.3-17                 sass_0.4.2                 
## [159] Rcpp_1.0.9                  bookdown_0.29              
## [161] mvtnorm_1.1-3               circlize_0.4.15            
## [163] FlowSOM_2.4.0               RProtoBufLib_2.8.0         
## [165] fansi_1.0.3                 tzdb_0.3.0                 
## [167] ModelMetrics_1.2.2.2        parallelly_1.32.1          
## [169] R6_2.5.1                    grid_4.2.1                 
## [171] ggridges_0.5.4              lifecycle_1.0.2            
## [173] zip_2.2.1                   curl_4.3.2                 
## [175] ggsignif_0.6.3              googlesheets4_1.0.1        
## [177] jquerylib_0.1.4             Matrix_1.5-1               
## [179] RcppAnnoy_0.0.19            TH.data_1.1-1              
## [181] RColorBrewer_1.1-3          iterators_1.0.14           
## [183] gower_1.0.0                 svgPanZoom_0.3.4           
## [185] htmlwidgets_1.5.4           beachmat_2.12.0            
## [187] polyclip_1.10-0             terra_1.6-17               
## [189] rvest_1.0.3                 ComplexHeatmap_2.12.1      
## [191] globals_0.16.1              codetools_0.2-18           
## [193] lubridate_1.8.0             randomForest_4.7-1.1       
## [195] metapod_1.4.0               gtools_3.9.3               
## [197] dbplyr_2.2.1                R.methodsS3_1.8.2          
## [199] gtable_0.3.1                DBI_1.1.3                  
## [201] httr_1.4.4                  highr_0.9                  
## [203] KernSmooth_2.23-20          stringi_1.7.8              
## [205] vroom_1.5.7                 reshape2_1.4.4             
## [207] farver_2.1.1                hexbin_1.28.2              
## [209] Rgraphviz_2.40.0            timeDate_4021.104          
## [211] magick_2.7.3                DT_0.25                    
## [213] xml2_1.3.3                  colorRamps_2.3.1           
## [215] ggcyto_1.24.1               BiocNeighbors_1.14.0       
## [217] interp_1.1-3                scattermore_0.8            
## [219] bit_4.0.4                   jpeg_0.1-9                 
## [221] ggraph_2.0.6                pkgconfig_2.0.3            
## [223] gargle_1.2.1                rstatix_0.7.0              
## [225] knitr_1.40
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