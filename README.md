[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6806449.svg)](https://doi.org/10.5281/zenodo.6806449)

# R based analysis workflow for IMC data

<!-- badges: start -->
[![build](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml/badge.svg)](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml)
<!-- badges: end -->

R workflow highlighting analyses approaches for imaging mass cytometry (IMC; or other multiplexed imaging) data.

## Scope


This workflow explains the use of common R/Bioconductor packages to pre-process and analyse single-cell data obtained from segmented IMC images.
While we use IMC data as an example, the concepts presented here can be applied to images obtained by other technologies (e.g. CODEX, MIBI, mIF, etc.).
The workflow can be largely divided into the following parts:

1. Preprocessing (reading in the data, spillover correction)
2. Image- and cell-level quality control, low-dimensional visualization
3. Sample/batch effect correction
4. Cell phenotyping via clustering or classification
5. Image visualization
6. Spatial analyses

## Usability

After cloning the repository, the code can be run as is.
It is continously tested on a ubuntu system using the newest release versions of the used R packages.

## Contribution

For feature requests and bug reports, please raise an issue [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/issues).
You can also add sections by forking the repo, adding your changes and issuing a pull request.

## Maintainer

[Nils Eling](https://github.com/nilseling)

## Contributors

[Vito Zanotelli](https://github.com/votti)
[Daniel Schulz](https://github.com/SchulzDan)
[Jonas Windhager](https://github.com/jwindhager) 
[Michelle Daniel](https://github.com/michdaniel)

## Citation

The workflow is currently under development for final publication.
In the meantime please refer to the 
[preprint](https://www.biorxiv.org/content/10.1101/2021.11.12.468357v1) 
which you can site as follows:

```
Jonas Windhager, Bernd Bodenmiller, Nils Eling (2020). An end-to-end workflow for multiplexed image processing and analysis. 
    bioRxiv, doi: 10.1101/2021.11.12.468357
```