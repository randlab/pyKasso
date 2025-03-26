![pyKasso's banner](/img/pykasso_banner_logo.png)

<!-- ![]() -->
<!-- [![PyPI Version](https://img.shields.io/pypi/v/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->
<!-- [![PyPI Status](https://img.shields.io/pypi/status/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->
<!-- [![PyPI Versions](https://img.shields.io/pypi/pyversions/pykasso.png)](https://pypi.python.org/pypi/pykasso) -->

![license](https://img.shields.io/github/license/randlab/pyKasso)
![last-commit](https://img.shields.io/github/last-commit/randlab/pyKasso/master)

<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/randlab/pyKasso/dev) -->

## pyKasso: a stochastic karst network simulation tool

pyKasso is a python3 open-source package intended to simulate easily and quickly karst networks using a geological model, hydrogeological, and structural data. It relies on a pseudo-genetic methodology where stochastic data and fast-marching methods are combined to perform thousands of simulations rapidly. The method is based on the stochastic karst simulator developed by Borghi et al (2012). It has been extended to account for anisotropy allowing to simplify the algorithm while accounting better for the geological structure following the method presented in Fandel et al. (2022). Statistical geometrical and topological metrics are computed on the simulated networks and compared with the same statistics computed on real karst network to evaluate the plausibility of the simulations.

![gif_01](/img/animation_01.gif)
![gif_02](/img/animation_02.gif)

## Installation

> [!IMPORTANT]
> Currently, pyKasso is working with:
> - Python 3.11 
> - Python 3.10 (not tested)
> - Python 3.9

The easiest way to install pyKasso is to use anaconda or miniconda.

1. Create a python 3.11 environnement:
```
conda create --name pyKasso -c conda-forge python=3.11
conda activate pyKasso
```

2. Clone the project then install it locally:
```
pip install -e git+https://github.com/randlab/pyKasso.git#egg=pykasso
```
Or download it then install it locally:
```
python -m pip install .
```

3. Install the hfm package:
```
conda config --add channels agd-lbr
conda install hfm
```

An alternative and faster way to do the same is to use miniforge: https://github.com/conda-forge/miniforge
```
conda create --name pyKasso -c conda-forge python=3.11
conda activate pyKasso
python -m pip install .
conda config --add channels agd-lbr
mamba install hfm
```

## Examples

- Examples developped for the paper: [notebooks/paper/](https://github.com/randlab/pyKasso/tree/master/notebooks/paper)
- Some basic examples illustrating pyKasso's functionalities: [notebooks/geometry/](https://github.com/randlab/pyKasso/tree/master/notebooks/geometry)
- An example to use pyKasso with Google Colab: [notebooks/colab/](https://github.com/randlab/pyKasso/tree/master/notebooks/colab) 

## Citing pyKasso

We published a [paper about pyKasso](https://doi.org/10.1016/j.envsoft.2025.106362).

If you are using pyKasso in your scientific research, please help our scientific visibility by citing our work.

> Miville, F., Renard, P., Fandel, C., Filipponi, M. 2025: pyKasso: An open-source three-dimensional discrete karst network generator. Environmental Modelling & Software, Volume 186, https://doi.org/10.1016/j.envsoft.2025.106362

BibTex:
```
@article{MIVILLE2025106362,
  title = {pyKasso: An open-source three-dimensional discrete karst network generator},
  journal = {Environmental Modelling & Software},
  volume = {186},
  pages = {106362},
  year = {2025},
  issn = {1364-8152},
  doi = {https://doi.org/10.1016/j.envsoft.2025.106362},
  url = {https://www.sciencedirect.com/science/article/pii/S1364815225000465},
  author = {François Miville and Philippe Renard and Chloé Fandel and Marco Filipponi},
}
```

## Publications

- Miville, F., Renard, P., Fandel, C., Filipponi, M. 2025: pyKasso: An open-source three-dimensional discrete karst network generator. Environmental Modelling & Software, Volume 186, https://doi.org/10.1016/j.envsoft.2025.106362
- Fandel, C., Miville, F., Ferré, T. et al. 2022: The stochastic simulation of karst conduit network structure using anisotropic fast marching, and its application to a geologically complex alpine karst system. Hydrogeol J 30, 927–946, https://doi.org/10.1007/s10040-022-02464-x
- Borghi, A., Renard, P., Jenni, S. 2012: A pseudo-genetic stochastic model to generate karstic networks, Journal of Hydrology, 414–415, https://doi.org/10.1016/j.jhydrol.2011.11.032.
