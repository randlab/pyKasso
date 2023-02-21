# pyKasso - Stochastic karst network simulation

pyKasso is a python3 project to stochastically simulate karst networks.

This is the stable 2D version as of June 2022.

## Installation
Note: This may not work on non-Windows machines

1. Clone this repository (make sure you are in the correct branch).
    - Click the green "Code" button on the top right
    - Clone using your method of choice (GitHub Desktop for example)
	- `git clone -b cfpy https://github.com/randlab/pyKasso.git`
2. Clone the karsnet repository and keep track of which directory it is in:
    - Karsnet GitHub page is [here](https://github.com/karstnet/karstnet)
	- `git clone https://github.com/karstnet/karstnet.git`
3. Create and activate a new environment called pykasso2D from the environment.yml file provided:
    - `conda env create --file your/path/to/environment.yml`
    - `conda activate pykasso2D`    
4. Manually install karstnet and pykasso from the folders you cloned:
    - `pip install -e your/path/to/karstnet`
	- `pip install -e your/path/to/pykasso`
5. Test by opening a Jupyter notebook, and running a simple example:
	- `jupyter notebook`
	- In Jupyter notebook, navigate in the directory to your GitHub folder and to the pyKasso repository. 
	- Find the folder called "notebooks"
	- Open the notebook "simple_fastmarching_example.ipynb"
	- Change the paths in the notebook to your local paths
	- Try to run the notebook

You may experience issues with the `agd-hfm` package. If such a message appear during pykasso execution
```
ERROR : the HFMpy library is not found, and the path to the FileHFM binaries is not set
```
it means you have to manually install the binaries from agd-hfm packgage, according to your python version and your operating system you will select the adequate binaries here : 
https://anaconda.org/AGD-LBR/hfm/files
You should unzip it twice, and directly add the content in the `site-packages` directory of your current python environnement.
You can easily locate the root of your environnement with `conda env list`, `site-packages` should be then in the `Lib/` directory.  


Note: You may need to install some of the packages listed in the dependencies below if they didn't get installed properly from the environment.yml file. 

OR, try the older simpler install version:

Install from source (you must first change directories to the project main directory):
```
pip install -e your/path/to/pykasso
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
- openpyxl
- [karstnet](https://github.com/UniNE-CHYN/karstnet) use `pip install -e your\path\to\karstnet`
	- networkx
	- scipy
	- mplstereonet (Note: cannot install with conda - use pip)

## Examples

The folder `notebooks/` has several examples:
- simple_fastmarching_example: shows the difference between isotropic and anisotropic fast marching on a square 9x9 grid
- simple_geology_generator: shows examples of different ways to load or randomly generate geologic settings for the karst network model (interacting with DEMs, GSLIB, VTK, etc.)
- Fandel_et_al_2022: shows the code used to generate figures for the paper Fandel et al. 2022, which demonstrates the capabilities of pyKasso in a simple valley, and then the application of pyKasso to a complex alpine karst catchment (Gottesacker)


## File structure

Examples of Jupyter Notebooks using pyKasso are in the notebooks folder. 
Within this folder, the inputs folder has a subfolder for each study site:
- Betteraz (Switzerland)
- Gottesacker (Germany/Austria)
- Tsanfleuron (Switzerland)
- Valley1 (simple synthetic example)