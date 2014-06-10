
[![Build Status](https://travis-ci.org/mwpennell/arbutus.png?branch=master)](https://travis-ci.org/mwpennell/arbutus)

![arbutus logo](https://github.com/mwpennell/arbutus/raw/master/extra/arbutus_logo.png)

This package was developed to assess the adequacy of continuous trait models. The package includes functions to parse the output of various model-fitting functions in R, rescale phylogenies according to fitted parameters, simulate datasets, calculate summary statistics and plot results. The methods is described in a forthcoming paper: "Model adequacy and the macroevolution of angiosperm functional traits" [(preprint here)](http://biorxiv.org/content/early/2014/04/07/004002).

## Installation

To install the package from github, the easiest way is to install using [devtools](https://github.com/hadley/devtools). Install `devtools`, then type

```
library(devtools)
install_github("mwpennell/arbutus")
```

The package depends on the `ape` library



## Supported packages

The idea behind `arbutus` is that models are fit to comparative data using other packages and `arbutus` is used to assess the adequacy of the fitted model. There are no native functions for fitting models or performing model selection within this package. Currently, arbutus can handle objects of the following classes:

* `gfit` objects from `fitContinuous` in [geiger](https://github.com/mwpennell/geiger-v2) (>= 2.0)
* `fit.mle` objects from `find.mle` in [diversitree](https://github.com/richfitz/diversitree) (>= 0.9.7)
* `mcmcsamples` objects from `mcmc` in [diversitree](https://githbu.com/richfitz/diversitree) (>= 0.9.7)`
* `gls` objects from `gls` in `nlme`
* `pgls` objects from `pgls` in `caper`
* `phylolm` objects from `phylolm` in `phylolm`
* `phylo` objects in `ape` format (does not rescale phylogeny)
* `multiPhylo` objects in `ape` format (does not rescale phylogeny)

We plan to add more package interfaces as this project progresses.

