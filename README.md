# R based analysis workflow for IMC data

<!-- badges: start -->
[![build](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml/badge.svg)](https://github.com/BodenmillerGroup/IMCDataAnalysis/actions/workflows/build.yml)
<!-- badges: end -->

R workflow highlighting analyses approaches for imaging mass cytometry (IMC; or other multiplexed imaging) data.

## Scope

**This project is under development and will change on a regular basis**

Upon release, this workflow highlights the use of common R/Bioconductor packages to pre-process and analyse single-cell data obtained from segmented IMC images.
While we use IMC data as an example, the concepts presented here can be applied to images obtained by other technologies (e.g. CODEX, MIBI, mIF, etc.).
The workflow can be largely divided into the following parts:

1. Preprocessing (reading in the data, spillover correctio)
2. Metrics for quality control
3. Low-dimensional visualization, clustering and/or cell-type classification
4. Visualization of cell- and pixel-level information
5. Spatial analyses

## Usability

After cloning the repository, the code can be run as is.
It is continously tested on a ubuntu system using the newest release versions of the used packages.

## Contribution

For feature requests and bug reports, please raise an issue [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/issues).
You can also add sections by forking the repo, adding your changes and issuing a pull request.

## Maintainer

[Nils Eling](https://github.com/nilseling)