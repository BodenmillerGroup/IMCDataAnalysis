[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6806449.svg)](https://doi.org/10.5281/zenodo.6806449)

# R based analysis workflow for multiplexed imaging data

<!-- badges: start -->
[![build](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml/badge.svg)](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml)
<!-- badges: end -->

R workflow highlighting analyses approaches for multiplexed imaging data.

## Scope


This workflow explains the use of common R/Bioconductor packages to pre-process and analyse single-cell data obtained from segmented multichannel images.
While we use imaging mass cytometry (IMC) data as an example, the concepts presented here can be applied to images obtained by other technologies (e.g. CODEX, MIBI, mIF, CyCIF, etc.).
The workflow can be largely divided into the following parts:

1. Preprocessing (reading in the data, spillover correction)
2. Image- and cell-level quality control, low-dimensional visualization
3. Sample/batch effect correction
4. Cell phenotyping via clustering or classification
5. Single-cell visualization
6. Image visualization
7. Spatial analyses

## Usage

To reproduce the analysis displayed at [https://bodenmillergroup.github.io/IMCDataAnalysis/](https://bodenmillergroup.github.io/IMCDataAnalysis/) clone the repository via:

```
git clone https://github.com/BodenmillerGroup/IMCDataAnalysis.git
```

For reproducibility purposes, we provide a Docker container [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/pkgs/container/imcdataanalysis).

1. After installing [Docker](https://docs.docker.com/get-docker/) you can run the container via:

```
docker run –v /path/to/IMCDataAnalysis:/home/rstudio/IMCDataAnalysis \
	-e PASSWORD=bioc -p 8787:8787  \
	ghcr.io/bodenmillergroup/imcdataanalysis:latest
```

2. An RStudio server session can be accessed via a browser at `localhost:8787` using `Username: rstudio` and `Password: bioc`.  
3. Navigate to `IMCDataAnalysis` and open the `IMCDataAnalysis.Rproj` file.  
4. Code in the individual files can now be executed or the whole workflow can be build by entering `bookdown::render_book()`.

## Contributing guidelines

For feature requests and bug reports, please raise an issue [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/issues).

For adding new content to the book please work inside the Docker container as explained above.
You can fork the repository, add your changes and open a pull request.
To add new libraries to the container please add them to the [Dockerfile](Dockerfile).

## Maintainer

[Nils Eling](https://github.com/nilseling)

## Contributors

[Vito Zanotelli](https://github.com/votti)  
[Daniel Schulz](https://github.com/SchulzDan)  
[Jonas Windhager](https://github.com/jwindhager)   
[Michelle Daniel](https://github.com/michdaniel)  
[Lasse Meyer](https://github.com/lassedochreden)

## Citation

The workflow is currently under development for final publication.
In the meantime please refer to the 
[preprint](https://www.biorxiv.org/content/10.1101/2021.11.12.468357v1) 
which you can site as follows:

```
Jonas Windhager, Bernd Bodenmiller, Nils Eling (2020). An end-to-end workflow for multiplexed image processing and analysis. 
    bioRxiv, doi: 10.1101/2021.11.12.468357
```
