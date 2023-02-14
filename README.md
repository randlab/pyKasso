![pyKasso's logo](/docs/source/_static/pykasso_logo.png)

## pyKasso: a stochastic karst network simulation tool

pyKasso is a python3 open-source package intended to simulate easily and quickly karst networks using a geological model, hydrogeological, and structural data. It relies on a pseudo-genetic methodology where stochastic data and fast-marching methods are combined to perform thousands of simulations rapidly. The method is based on the stochastic karst simulator developed by Borghi et al (2012). It has been extended to account for anisotropy allowing to simplify the algorithm while accounting better for the geological structure following the method presented in Fandel et al. (2022). Statistical geometrical and topological metrics are computed on the simulated networks and compared with the same statistics computed on real karst network to evaluate the plausibility of the simulations.

## Quick installation

Work in progress

### With pip

From source:
```
python -m pip install .
```

### With poetry

From source:
```
poetry install
```

### Check installation

```
poetry run pytest tests/
```

### Dependencies

pyKasso requires the following python packages to function properly:
- [agd](https://github.com/Mirebeau/AdaptiveGridDiscretizations)
- [karstnet](https://github.com/UniNE-CHYN/karstnet)

## Documentation

Work in progress

## Examples

Work in progress

## Contact

- F. Miville
- Prof. C. Fandel
- Prof. P. Renard

## Publications

- Fandel, C., Miville, F., Ferré, T. et al. 2022: The stochastic simulation of karst conduit network structure using anisotropic fast marching, and its application to a geologically complex alpine karst system. Hydrogeol J 30, 927–946, https://doi.org/10.1007/s10040-022-02464-x
- Borghi, A., Renard, P., Jenni, S. 2012: A pseudo-genetic stochastic model to generate karstic networks, Journal of Hydrology, 414–415, https://doi.org/10.1016/j.jhydrol.2011.11.032.