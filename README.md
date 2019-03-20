Modular R library for multiple Pathway Enrichment methods and their integration in KNIME using GenericKnimeNodes
================================================================================================================

The code can be cloned like this:
```
git clone https://github.com/InesAssum/EnrichmentNodes
```

Licensing:
---------

The EnrichmentNodes source code is published under a ...I have no idea... license. The licenses of the individual applications 
integrated into EnrichmentNodes can be found in LICENSE.

* KNIME: General Public License v3
* Docker: Apache License v2.0
* Generic KNIME Nodes: Apache License v2.0
* R/Rstudio: AGPL v3
* mono: https://github.com/mono/mono/blob/master/LICENSE
* MONA: academic, non-commercial use only


Prerequisites:
-------------

1) Docker (https://www.docker.com/products/docker-desktop)

and for using our nodes in KNIME:

2) KNIME Analytics platform, we suggest version 3.7 (https://www.knime.com/downloads)
3) place the `export_plugin/plugins/de.enrichment_1.2.0.0.jar` under
  - Windows: `C:/Programme/KNIME 3.7.0/plugins`
  - Mac: `/⁨Applications/⁨KNIME 3.7.0.app⁩/Contents/Eclipse⁩/plugins⁩`
  - Linux: ``
  

### for developers:

Apache Ant

KNIME SDK:
```
git clone https://github.com/knime/knime-sdk-setup
```

Generic Knime Nodes (GKN):
```
git clone https://github.com/genericworkflownodes/GenericKnimeNodes/
```


Really quick how-to:
-------------------

### Use our R framework with the enrich docker:
on the shell, change to your desired working directory and execute:
```
docker run -d -p 8787:8787 -v %cd%:/home/rstudio -e PASSWORD=password enrich:latest
```
visit `localholst:8787` to enjoy your interactive rstudio version!

### KNIME
* start KNIME 3.7 and enjoy our nodes under `Community Nodes/EnrichmentNodes`

### Develop your own nodes
* write an R script or commandline tool that extends our existing R functions
* include all necessary scripts/executables in the src/files/ folder
* modify the existing or create a new docker, make sure that the scripts are executable
* add a line at the end of the plugin.properties defining your tool and docker
* add a .ctd file for your node to the knime/descriptors/ folder
* build with 'ant' and import the plugin including your tool!
* detailed tutorials will come later!


### everything else: checkout
```
https://github.com/InesAssum/EnrichmentNodes/tutorials
```


Contents
--------

### src

All files (including Dockerfile) for our enrich docker image. The `files` folder contains all code for the modular R framework.


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
```
https://join.slack.com/t/knime-setup/shared_invite/enQtNTMwMzk4NDE0MjQ2LTMwZDhmYmI2YjdhNDdmNmFmNTY3MDE4NzZiOGY0MWJkNWE5OTkwNjRkYzI1YzYwYWM5YjRlZDZjODg0MmFmMzA
```



References:
-----------

* GSEA (Subramanian, Aravind, et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proceedings of the National Academy of Sciences 102.43 (2005): 15545-15550.)
* fgsea (Sergushichev, Alexey. "An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation." BioRxiv (2016): 060012.)
* MONA (Sass, Steffen, et al. "A modular framework for gene set analysis integrating multilevel omics data." Nucleic acids research 41.21 (2013): 9622-9633.)
* KNIME (Berthold, Michael R., et al. "KNIME-the Konstanz information miner: version 2.0 and beyond." AcM SIGKDD explorations Newsletter 11.1 (2009): 26-31.)
* GKN (https://github.com/genericworkflownodes/GenericKnimeNodes)


