# irtpar
Estimate Bayesian IRT Models in Parallel.

This packages allows you to IRT models on sub samples of the data, in parallel and then combine them as described in [Neiswanger, Wang, Xing (2014)](http://arxiv.org/abs/1311.4780).

## Installation

Some functionality can be used already. In order to be able to install the package you must have *rstan* installed. To do that follow the instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Then *irtpar* can be installed from within *R* with the following steps:
```{r}
library(devtools)
install_github("flinder/irtpar")
```

