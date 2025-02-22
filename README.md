# ![nfcore/imcyto](docs/images/nfcore-imcyto_logo.png)

[![Build Status](https://travis-ci.com/nf-core/imcyto.svg?branch=master)](https://travis-ci.com/nf-core/imcyto)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![Docker](https://img.shields.io/docker/automated/nfcore/imcyto.svg)](https://hub.docker.com/r/nfcore/imcyto)

## Introduction

**nfcore/imcyto** is a bioinformatics analysis pipeline used for Image Mass Cytometry analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline summary

1. Split mcd file by ROI, and save full and ilastik stacks separately based on specification in `metadata.csv` ([`imctools`](https://github.com/BodenmillerGroup/imctools))
2. Apply preprocessing filters to full stack tiff files ([`CellProfiler`](https://cellprofiler.org/); `--full_stack_cppipe` parameter)
3. Merge images from ilastik stack to obtain RGB image of cell nuclei and membranes to generate a composite tiff ([`CellProfiler`](https://cellprofiler.org/); `--ilastik_stack_cppipe` parameter)
4. Use composite tiff to classify pixels as membrane, nuclei or background, and save probabilities map as tiff ([`Ilastik`](https://www.ilastik.org/); *optional*)
5. Use probability tiffs and preprocessed full stack tiffs for single cell segmentation to generate a cell mask as tiff and then overlay cell mask onto full stack tiff images to extract single cell information generating a csv file ([`CellProfiler`](https://cellprofiler.org/); `--segmentation_cppipe` parameter)

## Documentation

The nf-core/imcyto pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits
The pipeline was originally written by the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) for use at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Harshil Patel](mailto:harshil.patel@crick.ac.uk) and [Nourdine Bah](mailto:nourdine.bah@crick.ac.uk) in collaboration with [Karishma Valand](mailto:karishma.valand@crick.ac.uk), [Febe van Maldegem](mailto:febe.vanmaldegem@crick.ac.uk), [Emma Colliver](mailto:emma.colliver@crick.ac.uk) and [Mihaela Angelova](mailto:mihaela.angelova@crick.ac.uk).

Many thanks to others who contributed as a result of the Crick Data Challenge (Jan 2019) - [Gavin Kelly](mailto:gavin.kelly@crick.ac.uk), [Becky Saunders](mailto:becky.saunders@crick.ac.uk), [Katey Enfield](mailto:katey.enfield@crick.ac.uk), [Alix Lemarois](mailto:alix.lemarois@crick.ac.uk), [Nuria Folguera Blasco](mailto:nuria.folguerablasco@crick.ac.uk), [Andre Altmann](mailto:a.altmann@ucl.ac.uk).

It would not have been possible to develop this pipeline without the guidelines, scripts and plugins provided by the [Bodenmiller Lab](http://www.bodenmillerlab.com/). Thank you too!

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/imcyto for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
