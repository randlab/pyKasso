![pyKasso's banner](/docs/source/_static/pykasso_banner_logo.png)

<!-- ![]() -->
[![PyPI Version](https://img.shields.io/pypi/v/pykasso.png)](https://pypi.python.org/pypi/pykasso)
[![PyPI Status](https://img.shields.io/pypi/status/pykasso.png)](https://pypi.python.org/pypi/pykasso)
[![PyPI Versions](https://img.shields.io/pypi/pyversions/pykasso.png)](https://pypi.python.org/pypi/pykasso)

![license](https://img.shields.io/github/license/randlab/pyKasso)
![last-commit](https://img.shields.io/github/last-commit/randlab/pyKasso/dev)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/randlab/pyKasso/dev)

## pyKasso: a stochastic karst network simulation tool
<!-- ![pyKasso's logo](/docs/source/_static/pykasso_logo.png) -->

pyKasso is a python3 open-source package intended to simulate easily and quickly karst networks using a geological model, hydrogeological, and structural data. It relies on a pseudo-genetic methodology where stochastic data and fast-marching methods are combined to perform thousands of simulations rapidly. The method is based on the stochastic karst simulator developed by Borghi et al (2012). It has been extended to account for anisotropy allowing to simplify the algorithm while accounting better for the geological structure following the method presented in Fandel et al. (2022). Statistical geometrical and topological metrics are computed on the simulated networks and compared with the same statistics computed on real karst network to evaluate the plausibility of the simulations.

![gif_01](/docs/source/_static/animation_01.gif)
![gif_02](/docs/source/_static/animation_02.gif)

## Installation

Currently, pyKasso is only working with Python 3.9.

### Using conda

Download *environment.yml*. From source:
```
conda env create --name pykasso --file=environment.yml
```

<!-- Then:
```
pip install -e pykasso[analysis, visualization]
``` -->

### Using pip

In local, from source:
```
pip install -e ".[analysis,visualization,trame]"
```

With git:
```
pip install -e git+https://github.com/randlab/pyKasso.git@dev#egg=pykasso[analysis,visualization,trame]
```

<!-- 
### Check installation

Work in progress.

```
poetry run pytest tests/
```


### Dependencies

pyKasso requires the following python packages to function properly:
- [agd](https://github.com/Mirebeau/AdaptiveGridDiscretizations)
- [karstnet](https://github.com/UniNE-CHYN/karstnet)
- [pyvista](https://github.com/pyvista/pyvista)
-->

## Documentation

Work in progress.

## Examples

Some basic examples are avaible here : [notebooks/geometry/](https://github.com/randlab/pyKasso/tree/dev/notebooks/geometry)
Other examps : [notebooks/paper/](https://github.com/randlab/pyKasso/tree/dev/notebooks/paper)

## Contact

- F. Miville
- Prof. C. Fandel
- Prof. P. Renard

## Publications

- Fandel, C., Miville, F., Ferré, T. et al. 2022: The stochastic simulation of karst conduit network structure using anisotropic fast marching, and its application to a geologically complex alpine karst system. Hydrogeol J 30, 927–946, https://doi.org/10.1007/s10040-022-02464-x
- Borghi, A., Renard, P., Jenni, S. 2012: A pseudo-genetic stochastic model to generate karstic networks, Journal of Hydrology, 414–415, https://doi.org/10.1016/j.jhydrol.2011.11.032.