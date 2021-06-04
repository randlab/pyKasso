.. pyKasso documentation master file, created by
   sphinx-quickstart on Mon May 31 11:32:28 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyKasso's documentation
***********************

pyKasso is a python 3 package intended to simulate stochastic karst network.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Installation
============

pyKasso requires several other packages to properly works. 
The easiest way to install pykasso is using Anaconda.

In the Anaconda prompt:
```
conda env update -n pykasso --file environment.yaml
```

Then install pykasso from source (from project main directory):
```
pip install .
```

Dependencies
------------

pyKasso requires the following python packages to function properly:
- yaml
- mpmath
- numpy
- pandas
- matplotlib
- [skfmm](https://github.com/scikit-fmm/scikit-fmm)
- [karstnet](https://github.com/UniNE-CHYN/karstnet)

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
- [karstnet](https://github.com/UniNE-CHYN/karstnet) use `pip install -e your\path\to\karstnet`
	- networkx
	- scipy
	- mplstereonet (Note: cannot install with conda - use pip)
	
Examples
========

In ```notebooks```, several files have been provided in order to illustrate :
- How the classes work and interact between them;
- How to set up a model and start karst network modeling; 

Authors
=======

- F. Miville
- C. Fandel
- P. Renard

References
==========

- Article de Chloé 
- Article d'Andreas
- https://github.com/Mirebeau/HamiltonFastMarching

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
