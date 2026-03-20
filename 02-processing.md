# Multi-channel image processing {#processing}

This book focuses on common analysis steps of spatially-resolved single-cell data
**after** image segmentation and feature extraction. In this chapter, the sections
describe the processing of multiplexed imaging data, including file type
conversion, image segmentation, feature extraction and data export. To obtain
more detailed information on the individual image processing approaches, please
visit their repositories:

[steinbock](https://github.com/BodenmillerGroup/steinbock): The `steinbock`
toolkit offers tools for multi-channel image processing using the command-line
or Python code [@Windhager2021]. Supported tasks include IMC data pre-processing,
multi-channel image segmentation, object quantification and data
export to a variety of file formats. It supports functionality similar to those
of the IMC Segmentation Pipeline (see below) and further allows deep-learning enabled image
segmentation. The toolkit is available as platform-independent Docker
container, ensuring reproducibility and user-friendly installation. Read more in
the [Docs](https://bodenmillergroup.github.io/steinbock/latest/).

[IMC Segmentation
Pipeline](https://github.com/BodenmillerGroup/ImcSegmentationPipeline): The IMC
segmentation pipeline offers a rather manual way of segmenting multi-channel
images using a pixel classification-based approach. We continue to maintain the
pipeline but recommend the use of the `steinbock` toolkit for multi-channel
image processing.  Raw IMC data pre-processing is performed using the
[readimc](https://github.com/BodenmillerGroup/readimc) Python package to convert
raw MCD files into OME-TIFF and TIFF files. After image cropping, an
[Ilastik](https://www.ilastik.org/) pixel classifier is trained for image
classification prior to image segmentation using
[CellProfiler](https://cellprofiler.org/). Features (i.e., mean pixel intensity)
of segmented objects (i.e., cells) are quantified and exported. Read more in the
[Docs](https://bodenmillergroup.github.io/ImcSegmentationPipeline/).

## Image pre-processing (IMC specific)

Image pre-processing is technology dependent. While most multiplexed imaging
technologies generated TIFF or OME-TIFF files which can be directly segmented
using the `steinbock` toolkit, IMC produces data in the proprietary
data format MCD. 

To facilitate IMC data pre-processing, the
[readimc](https://github.com/BodenmillerGroup/readimc) open-source Python
package allows extracting the multi-modal (IMC acquisitions, panoramas),
multi-region, multi-channel information contained in raw IMC images. Both the
IMC Segmentation Pipeline and the `steinbock` toolkit use the `readimc`
package for IMC data pre-processing. Starting from IMC raw data and a "panel"
file, individual acquisitions are extracted as TIFF files and OME-TIFF files if
using the IMC Segmentation Pipeline. The panel contains information of
antibodies used in the experiment and the user can specify which channels to
keep for downstream analysis. When using the IMC Segmentation Pipeline, random
tiles are cropped from images for convenience of pixel labelling.

## Image segmentation

The IMC Segmentation Pipeline supports pixel classification-based image
segmentation while `steinbock` supports pixel classification-based and deep
learning-based segmentation.

**Pixel classification-based** image segmentation is performed by training a 
random forest classifier using [Ilastik](https://www.ilastik.org/) on the
randomly extracted image crops and selected image channels. Pixels are
classified as nuclear, cytoplasmic, or background. Employing a customizable
[CellProfiler](https://cellprofiler.org/) pipeline, the probabilities are then
thresholded for segmenting nuclei, and nuclei are expanded into cytoplasmic
regions to obtain cell masks.

**Deep learning-based** image segmentation is performed as presented by
[@Greenwald2021]. Briefly, `steinbock` first aggregates user-defined
image channels to generate two-channel images representing nuclear and
cytoplasmic signals. Next, the
[DeepCell](https://github.com/vanvalenlab/intro-to-deepcell) Python package is
used to run `Mesmer`, a deep learning-enabled segmentation algorithm pre-trained
on `TissueNet`, to automatically obtain cell masks without any further user
input.

Segmentation masks are single-channel images that match the input images in
size, with non-zero grayscale values indicating the IDs of segmented objects
(e.g., cells). These masks are written out as TIFF files after segmentation.

## Feature extraction {#feature-extraction}

Using the segmentation masks together with their corresponding multi-channel
images, the IMC Segmentation Pipeline as well as the `steinbock` toolkit extract
object-specific features. These include the mean pixel intensity per object and
channel, morphological features (e.g., object area) and the objects' locations.
Object-specific features are written out as CSV files where rows represent
individual objects and columns represent features.

Furthermore, the IMC Segmentation Pipeline and the `steinbock` toolkit compute
_spatial object graphs_, in which nodes correspond to objects, and nodes in
spatial proximity are connected by an edge. These graphs serve as a proxy for
interactions between neighboring cells. They are stored as edge list in form of
one CSV file per image.

Both approaches also write out image-specific metadata (e.g., width and height)
as a CSV file.

## Data export

To further facilitate compatibility with downstream analysis, `steinbock`
exports data to a variety of file formats such as OME-TIFF for images, FCS for
single-cell data, the _anndata_ format [@Virshup2021] for data analysis in Python,
and various graph file formats for network analysis using software such as
[CytoScape](https://cytoscape.org/) [@Shannon2003]. For export to OME-TIFF,
steinbock uses [xtiff](https://github.com/BodenmillerGroup/xtiff), a Python
package developed for writing multi-channel TIFF stacks.

## Data import into R

In Section \@ref(read-data), we will highlight the use of the
[imcRtools](https://github.com/BodenmillerGroup/imcRtools) and
[cytomapper](https://github.com/BodenmillerGroup/cytomapper) R/Bioconductor
packages to read spatially-resolved, single-cell and images as generated by the
IMC Segmentation Pipeline and the `steinbock` toolkit into the statistical
programming language R. All further downstream analyses are performed in R and
detailed in the following sections.






