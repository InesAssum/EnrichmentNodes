# Modular R library for multiple Pathway Enrichment methods and their integration in KNIME using GenericKnimeNodes


The code can be cloned like this:
```
git clone https://github.com/InesAssum/EnrichmentNodes
```


References: (TODO)
* GSEA
* fgsea
* MONA
* KNIME
* GKN
* Benjamin Schubert: Immunonodes

## Contents

### prerequisites

#### the KNIME SDK for developers
```
git clone https://github.com/knime/knime-sdk-setup
```

#### our nodes are created using Generic Knime Nodes (GKN)
```
git clone https://github.com/genericworkflownodes/GenericKnimeNodes/
```

#### for a nice example, also checkout Benjamin Schubert's ImmunoNodes:
```
git clone https://github.com/FRED-2/ImmunoNodes/
```


### src

All the code for the modular R framework can be found here. This also includes the Dockerfile used to setup our "enrich" container, which will later also be available at dockerhub.io.


### knime

This directory contains all the files needed for building the EnrichmentNodes via ant.

#### descriptors:
contain one .ctd file per node that determines the graphical interface and the input/output ports

#### plugin.properties
define the details of your plugin and the dockers used for executing each tool

#### DESCRIPTION, COPYRIGHT, LICENSE

#### You want to include your own tool?
* write an R script or commandline tool that extends our existing R functions
* include all necessary scripts/executables in the src/files/ folder
* modify the existing or create a new docker, make sure that the scripts are executable
* add a line at the end of the plugin.properties defining your tool and docker
* add a .ctd file for your node to the knime/descriptors/ folder
* build with 'ant' and import the plugin including your tool!
* detailed tutorials will come later!


### examples

Here you can find example data and script to toy around with in R as well as workflows you can try in KNIME.


### tutorials

Goal of the project is to create really detailed tutorials starting from scratch for
* methods / best practices / parameter description
* KNIME setup and node development
* hopefully all the stuff, that is totally logical now but took us weeks to figure out!


## You need help?

Come join our slack workspace!
```
https://join.slack.com/t/knime-setup/shared_invite/enQtNTMwMzk4NDE0MjQ2LTMwZDhmYmI2YjdhNDdmNmFmNTY3MDE4NzZiOGY0MWJkNWE5OTkwNjRkYzI1YzYwYWM5YjRlZDZjODg0MmFmMzA
```




