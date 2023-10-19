[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8100220.svg)](https://doi.org/10.5281/zenodo.6806448)

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

**Of note: it is recommended to use a date-tagged version of the container to ensure reproducibility**. 
This can be done via:

```
docker pull ghcr.io/bodenmillergroup/imcdataanalysis:<year-month-date>
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

Please cite the following paper when using the presented workflow in your research:

>  Windhager, J., Zanotelli, V.R.T., Schulz, D. et al. An end-to-end workflow for multiplexed image processing and analysis. Nat Protoc (2023). https://doi.org/10.1038/s41596-023-00881-0

    @article{Windhager2023,
        author = {Windhager, Jonas and Zanotelli, Vito R.T. and Schulz, Daniel and Meyer, Lasse and Daniel, Michelle and Bodenmiller, Bernd and Eling, Nils},
        title = {An end-to-end workflow for multiplexed image processing and analysis},
        year = {2023},
        doi = {10.1038/s41596-023-00881-0},
        URL = {https://www.nature.com/articles/s41596-023-00881-0},
        journal = {Nature Protocols}
    }


## Funding

The work was funded by the European Union’s Horizon 2020 research and innovation program under Marie Sklodowska-Curie Actions grant agreement No 892225 (N.E) and by the CRUK IMAXT Grand Challenge (J.W.).
