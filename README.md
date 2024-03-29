Modular R library for multiple Pathway Enrichment methods and their integration in KNIME using GenericKnimeNodes
================================================================================================================

The code can be cloned like this:
```
git clone https://github.com/InesAssum/EnrichmentNodes
```

Licensing:
---------

The EnrichmentNodes source code is published under a GPL3 license. The licenses of the individual applications 
integrated into EnrichmentNodes can be found in LICENSE and below.

* KNIME: General Public License v3
* Docker: Apache License v2.0
* Generic KNIME Nodes: Apache License v2.0
* R/Rstudio: AGPL v3
* mono: https://github.com/mono/mono/blob/master/LICENSE
* MONA: academic, non-commercial use only (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3834824/)
Sass S, Buettner F, Mueller NS, Theis FJ. A modular framework for gene set analysis integrating multilevel omics data. Nucleic Acids Res. 2013 Nov;41(21):9622-33. doi: 10.1093/nar/gkt752. Epub 2013 Aug 23. PubMed PMID: 23975194; PubMed Central PMCID: PMC3834824.


Prerequisites:
-------------

1) Docker (https://www.docker.com/products/docker-desktop)

and for using our nodes in KNIME:

2) KNIME Analytics platform, we suggest version 4.3.4 or higher (https://www.knime.com/downloads)
3) place the `export_plugin/en.enrichment_2.0.1.2.jar` under
  - Windows: `C:/Programme/KNIME 4.3.4/plugins`
  - Mac: `/⁨Applications/⁨KNIME 4.3.4.app⁩/Contents/Eclipse⁩/plugins⁩`
  - Linux: ``
  

### for developers:

Apache Ant

KNIME SDK:
```
https://github.com/knime/knime-sdk-setup
```
Import the KNIME SDK via Eclipse, the current nodes were generated using the branch `knime-sdk-setup releases/2020-12` and the `KNIME-AP-complete.target` platform run together with the `JavaSE-1.8 (AdoptOpenJDK (OpenJ9) 8 [1.8.0_275])` environment.

Generic Knime Nodes (GKN):
```
git clone https://github.com/genericworkflownodes/GenericKnimeNodes/
```


Really quick how-to:
-------------------

### Use our R framework with the enrich docker:
on the shell, change to your desired working directory and execute:

on mac/linux:
```
docker run -d -p 8787:8787 -v $(pwd):/home/rstudio -e PASSWORD=password inesassum/enrich:latest
```
or on Windows
```
docker run -d -p 8787:8787 -v %cd%:/home/rstudio -e PASSWORD=password inesassum/enrich:latest
```
visit `localhost:8787` in your browser to enjoy your interactive rstudio session with the `user: rstudio` and `pw: password`!

Older versions are available on docker.hub (https://hub.docker.com/repository/docker/inesassum/enrich)

### KNIME
* start KNIME 4.3.4 or later, drop our plugin in your plugin folder and enjoy our nodes under `Community Nodes/EnrichmentNodes`

### Develop your own nodes
* write an R script or commandline tool that extends our existing R functions
* include all necessary scripts/executables in the src/files/ folder
* modify the existing or create a new docker, make sure that the scripts are executable
* add a line at the end of the plugin.properties defining your tool and docker
* add a .ctd file for your node to the knime/descriptors/ folder
* build with 'ant' and import the plugin including your tool!
* detailed tutorials will come later!


### everything else: checkout
[Our tutorials](tutorials/)


Contents
--------

### src

All files (including Dockerfile) for our enrich docker image. The `library` folder contains all code for the modular R framework.


### knime

This directory contains all the files needed for building the EnrichmentNodes via ant.

#### descriptors:
contain one .ctd file per node that determines the graphical interface and the input/output ports

#### plugin.properties
define the details of your plugin and the dockers used for executing each tool

#### DESCRIPTION

#### COPYRIGHT

#### LICENSE


### examples

Here you can find example data and scripts to toy around with in R as well as workflows you can try in KNIME.


### tutorials

Goal of the project is to create really detailed tutorials starting from scratch for
* methods / best practices / parameter description
* KNIME setup and node development
* hopefully all the stuff, that is totally logical now but took us weeks to figure out!

for a nice example, also checkout Benjamin Schubert's ImmunoNodes:
```
git clone https://github.com/FRED-2/ImmunoNodes/
```

Need help?
----------

Come join our slack workspace!

https://join.slack.com/t/knime-setup/shared_invite/enQtNTMwMzk4NDE0MjQ2LTMwZDhmYmI2YjdhNDdmNmFmNTY3MDE4NzZiOGY0MWJkNWE5OTkwNjRkYzI1YzYwYWM5YjRlZDZjODg0MmFmMzA




References:
-----------

* GSEA (Subramanian, Aravind et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proceedings of the National Academy of Sciences 102.43 (2005): 15545-15550.)
* fgsea (Korotkevich, Gennady et al. (2021). "Fast gene set enrichment analysis." bioRxiv. doi: 10.1101/060012, https://www.biorxiv.org/content/10.1101/060012v3)
* MONA (Sass, Steffen et al. "A modular framework for gene set analysis integrating multilevel omics data." Nucleic acids research 41.21 (2013): 9622-9633.)
* MGSA (Bauer, Sebastian et al. "GOing Bayesian: model-based gene set analysis of genome-scale data." Nucleic acids research vol. 38,11 (2010): 3523-32. doi:10.1093/nar/gkq045)
* KNIME (Berthold, Michael R., et al. "KNIME-the Konstanz information miner: version 2.0 and beyond." AcM SIGKDD explorations Newsletter 11.1 (2009): 26-31.)
* GKN (https://github.com/genericworkflownodes/GenericKnimeNodes)


