# pyKasso - stochastic karst network simulation

pyKasso is a python3 project intended to simulate stochastic karst network.

Version 0.2.0 - May 2021

## Authors

François Miville <francois.miville@unine.ch> - Code architecture and data model
Chloé Fandel     <cfandel@email.arizona.edu> - Anisotropic FMM implementation
Philippe Renard  <philippe.renard@unine.ch>  - Fracture discretizer

## Installation
You can install the requirements for pykasso with the environment.yaml file:
```
conda env update -n pykasso --file environment.yaml
```

Then install pykasso from source (from project main directory):
```
pip install .
```

### Dependencies

pyKasso requires the following python packages to function properly:
- [agd-hfm](https://github.com/Mirebeau/AdaptiveGridDiscretizations)
	- On non-windows machines or for python 3.7, use `conda install agd -c agd-lbr`
	- On windows with python 3.8+:
		- clone the agd repository from GitHub
		- create a dedicated environment using the agd-hfm.yaml file:
		`conda env create --file path/to/agd.hfm.yaml`
		`conda activate agd-hfm.yaml`
		- install all other packages into this environment
		- note: you can change the name of the environment by changing the name of the agd-hfm.yaml file
- yaml (Note: use `conda install pyyaml` (not yaml))
- mpmath
- numpy
- pandas
- matplotlib
- xlrd
- [skfmm](https://github.com/scikit-fmm/scikit-fmm) use `conda install scikit-fmm`
- [karstnet](https://github.com/UniNE-CHYN/karstnet) use `pip install -e your\path\to\karstnet`
	- networkx
	- scipy
	- mplstereonet (Note: cannot install with conda - use pip)

To be able to use pykasso with swmmpy, the following packages are also needed:
- rasterio (Note: this will also install gdal? But may need to do a separate `conda install gdal` first)
- shapely
- pyproj
- pysheds
- swmmtoolbox (Note: cannot install with conda - use pip)




## Examples

Please have a look on the notebooks examples in ``notebooks/``.
