# pyKasso - stochastic karst network simulation

pyKasso is a python3 project intended to simulate stochastic karst network.

Version 0.1.0 - April 2020


## Installation

Install from source (from project main directory):
```
pip install .
```

# Dependencies

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


## Examples

Please have a look on the notebooks examples in ``notebooks/``.