# Useful information when developing this book

This document is to guide future developers to maintain and extend the IMC
data analysis book. 

## General setup

* The IMC data analysis book is written in [bookdown](https://bookdown.org/). 
* Each section is stored in its own `.Rmd` file with `index.Rmd` building the landing page
* References are stored in `book.bib`
* At the end of each `.Rmd` file a number of unit tests are executed. These 
unit tests are always executed but their results are not shown in the book.

### Continous integration/continous deployment

* CI/CD is executed based on the workflow [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/blob/main/.github/workflows/build.yml).
* On the first of each month based on the [Dockerfile](https://github.com/BodenmillerGroup/IMCDataAnalysis/blob/main/Dockerfile) a new Docker image is build. We are doing this so that the workflow is always tested against the newest software versions.
* The Docker image is pushed to the Github Container Registry [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/pkgs/container/imcdataanalysis).
* The Docker image is date tagged and `latest` always refers to the newest build.
* Once the Docker image is build, the IMC data analysis book is executed within the
newest Docker image. This will also run all unit tests.

**Of note:** Sometimes the calculation of the UMAP produces slightly different
results. If that happens the workflow run can be re-executed by clicking the `Re-run jobs` button of the workflow run.
This test could also be excluded on the long run.

* When pushing to `main` (either directly or via a PR), the CI/CD workflow is
executed. 
* If the Dockerfile changed (e.g., if you want to add a new package), a new Docker image is build and the workflow is executed within the new Docker image.
* If the Dockerfile did not change, the workflow is executed within the most recent Docker image.

## Updating the book

This section describes how to update the book. You want to do this to add new content
but also to fix bugs or adjust unit tests.

### Work on the devel branch

It is recommended to work on the `devel` branch of the Github repository to add
new changes. 

### Work within the newest Docker container

It is also recommended to always work within a Docker container based on the newest
Docker image available:

1. After installing [Docker](https://docs.docker.com/get-docker/) you can first pull the container via:

```
docker pull ghcr.io/bodenmillergroup/imcdataanalysis:yyyy-mm-dd
```

and then run the container:

```
docker run -v /path/to/IMCDataAnalysis:/home/rstudio/IMCDataAnalysis \
	-e PASSWORD=bioc -p 8787:8787  \
	ghcr.io/bodenmillergroup/imcdataanalysis:yyyy-mm-dd
```

2. An RStudio server session can be accessed via a browser at `localhost:8787` using `Username: rstudio` and `Password: bioc`.  
3. Navigate to `IMCDataAnalysis` and open the `IMCDataAnalysis.Rproj` file.  
4. Code in the individual files can now be executed or the whole workflow can be build by entering `bookdown::render_book()`.

### Adding new packages

If you need to add new packages to the workflow, make sure to add them to the
[software requirements](https://bodenmillergroup.github.io/IMCDataAnalysis/prerequisites.html#software-requirements)
section and to the Dockerfile.

### Opening a pull request

Now you can change the content of the book. 
Once you have added all changes, push the changes to `devel` and open a pull request
to `main`. Wait until all checks have passed and you can merge the PR.

### Add changes to CHANGELOG.md

Please track the changes that you are making in the [CHANGELOG.md](CHANGELOG.md) file.

### Trigger a new release

Once you have added the changes to the CHANGELOG, merged the pull request and 
the workflow has been executed on CI/CD, you can trigger a new release.

* Go to [here](https://github.com/BodenmillerGroup/IMCDataAnalysis/releases) and click on `Draft a new release` at the top of the page.
* Under `Choose a tag` create a new tag and give details on the release.
* With each release the corresponding [Zenodo repository](https://zenodo.org/records/10209942) is updated.

## Updating the data

For new `steinbock` releases and specifically if the Mesmer version changes, the
example data should be updated. The example data are stored on Central NAS 
and are hosted on Zenodo. 

### Re-analyse the example data

* You can find the raw data on [zenodo](https://zenodo.org/records/7575859).
* On Central NAS under projects/IMCWorkflow/zenodo create a new folder called `steinbock_0.x.y` where x denotes the new major version and y the new minor version.
* Copy the `steinbock.sh` script from the folder of the previous version to to folder of the newest version.
* Change the steinbock version number in the `steinbock.sh` script and execute it.
* It should generate all relevant files and zip all folders.

### Upload data to zenodo

* On [zenodo](https://zenodo.org/records/7624451), click on `New version` and replace all files with the newer version. No need to upload the raw data to zenodo as they are hosted in a different repository

### Adjust the book

* Work in the most recent Docker container and on the devel branch.
* Manually go through each section, update the links in the [Prerequisites](https://bodenmillergroup.github.io/IMCDataAnalysis/prerequisites.html#download-data) section
* Make sure to check and asjust the unit tests at the end of each file
* Make sure that the text (e.g. clustering) still matches the results

*Important:* as we are training a random forest classifier on manually gated cells, these gated cells won't match the newest version of the data if the Mesmer version changed. For this, we have the  `code/transfer_labels.R` script that automatically re-gates cells in the new SPE object.

* Go through all sections until `Cell phenotyping`
* Based on the old `gated_cells` and the new SPE object, execute the `code/transfer_labels.R` script
* Zip the new `gated_cells` and upload them to a new version on [zendod](https://zenodo.org/records/8095133)
* Adjust the link to the new gated cells in the [Prerequisites](https://bodenmillergroup.github.io/IMCDataAnalysis/prerequisites.html#download-data) section
* Make sure that the new classification results closely match the new results

* Continue going through the book

### Add changes to CHANGELOG.md

Finally, add all the recent changes to the CHANGELOG, create and merge a PR and create a new release (see above).


