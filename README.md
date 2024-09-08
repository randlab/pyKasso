![pyKasso's banner](/img/pykasso_banner_logo.png)

<!-- ![]() -->
<!-- [![PyPI Version](https://img.shields.io/pypi/v/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->
<!-- [![PyPI Status](https://img.shields.io/pypi/status/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->
<!-- [![PyPI Versions](https://img.shields.io/pypi/pyversions/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->

![license](https://img.shields.io/github/license/randlab/pyKasso)
![last-commit](https://img.shields.io/github/last-commit/randlab/pyKasso/dev)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/randlab/pyKasso/dev)

## pyKasso: a stochastic karst network simulation tool

pyKasso is a python3 open-source package intended to simulate easily and quickly karst networks using a geological model, hydrogeological, and structural data. It relies on a pseudo-genetic methodology where stochastic data and fast-marching methods are combined to perform thousands of simulations rapidly. The method is based on the stochastic karst simulator developed by Borghi et al (2012). It has been extended to account for anisotropy allowing to simplify the algorithm while accounting better for the geological structure following the method presented in Fandel et al. (2022). Statistical geometrical and topological metrics are computed on the simulated networks and compared with the same statistics computed on real karst network to evaluate the plausibility of the simulations.

![gif_01](/img/animation_01.gif)
![gif_02](/img/animation_02.gif)

## Installation

Currently, pyKasso is only working with Python 3.9. 

If you are using anaconda or miniconda, the simplest way to install pykasso is to create the `pykasso` conda environment av
nd activate it with the following commands:
```
conda env create -f environment.yml
conda activate pykasso
```
This will install all the required dependencies.

An alternative and faster way to do the same is to use miniforge: https://github.com/conda-forge/miniforge:
```
mamba env create -f environment.yml
mamba activate pykasso
```


## Examples

- Examples developped for the paper: [notebooks/paper/](https://github.com/randlab/pyKasso/tree/dev/notebooks/paper)
- Some basic examples illustrating pyKasso's functionalities: [notebooks/geometry/](https://github.com/randlab/pyKasso/tree/dev/notebooks/geometry)

## Publications

- Fandel, C., Miville, F., Ferré, T. et al. 2022: The stochastic simulation of karst conduit network structure using anisotropic fast marching, and its application to a geologically complex alpine karst system. Hydrogeol J 30, 927–946, https://doi.org/10.1007/s10040-022-02464-x
- Borghi, A., Renard, P., Jenni, S. 2012: A pseudo-genetic stochastic model to generate karstic networks, Journal of Hydrology, 414–415, https://doi.org/10.1016/j.jhydrol.2011.11.032.
