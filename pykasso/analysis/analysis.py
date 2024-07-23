"""
This module contains a class able to mangage project results in order to
perform statistical analysis.
"""

### Internal dependencies
import copy

### External dependencies
import numpy as np
import pandas as pd

### Optional dependencies
try:
    import karstnet as kn
except ImportError:
    _has_karstnet = False
else:
    _has_karstnet = True

### Typing
from typing import Union
from pykasso.core.project import Project
from pandas import (DataFrame, Series)
from pandas.io.formats.style import Styler


def requires_karstnet():
    """
    If 'karstnet' package is not installed, return `ImportError` exception
    when a method requiring 'karstnet' is called.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            if not _has_karstnet:
                msg = ("karstnet package is required to do this."
                       " 'pip install -e pykasso[analysis]' to install it.")
                raise ImportError(msg)
            result = function(*args, **kwargs)
            return result
        return _wrapper
    return _


class Analyzer():
    """
    This class manages pyKasso's project and provides methods to compute
    statistical analysis.
    """
    def __init__(self,
                 project: Project
                 ) -> None:
        """
        Initialize the class.
        """
        
        # Intialization
        self.project = project
        self.stats = None
        
        # Load reference metrics for statistical karstic network analysis
        self._load_statistics()
        
    def _load_statistics(self) -> None:
        """
        Set the reference metrics for statistical karstic network analysis
        More details here : https://github.com/karstnet/karstnet
        """
        package_location = self.project._pckg_paths['package_location']
        statistics_file_path = "/../_misc/statistics.xlsx"
        statistics_file_location = package_location + statistics_file_path
        self.stats = pd.read_excel(statistics_file_location).describe()
        return None
    
    @requires_karstnet()
    def compute_metrics(self,
                        verbose: bool = False,
                        ) -> DataFrame:
        """
        Compute the statistical metrics for each simulated karst network
        using the karstnet package.
        
        Parameters
        ----------
        verbosity : int, optional
            Verbosity of karstnet results, by default 0

        Returns
        -------
        df_metrics : pandas.DataFrame
            Dataframe of karstnet metrics.

        Notes
        -----
        Karstnet is a python3 project providing tools for the statistical
        analysis of karstic networks. More details here:
        https://github.com/karstnet/karstnet

        References
        ----------
        .. [1] Collon, P., Bernasconi D., Vuilleumier C., and Renard P., 2017,
               Statistical metrics for the characterization of karst network
               geometry and topology. Geomorphology. 283: 122-142 doi:10.1016/
               j.geomorph.2017.01.034
               http://dx.doi.org/doi:10.1016/j.geomorph.2017.01.034

        .. warning::
            A corrigendum has been published in Geomorphology journal:
            Geomorphology 389, 107848,
            http://dx.doi.org/doi:10.1016/j.geomorph.2021.107848.

        Examples
        --------
        >>> app = pk.pykasso()
        >>> TODO
        >>> df_metrics = app.analyzer.compute_metrics()
        """
        df_metrics = pd.DataFrame()

        # For each simulation, retrieve data and compute metrics
        for i, data in enumerate(self.project):
            
            # Retrieve data
            karstnet_edges = data["vectors"]["edges_"].to_numpy().tolist()
            karstnet_nodes = copy.deepcopy(data["vectors"]["nodes_"])
            
            # Drop last item in list (the node type) for each dictionary entry
            karstnet_nodes = karstnet_nodes.drop(columns=['type', 'vadose'])
            index = karstnet_nodes.index
            values = karstnet_nodes.iloc[:, [0, 1, 2]].to_numpy().tolist()
            karstnet_nodes = {i: value for i, value in zip(index, values)}

            # Compute karstnet metrics
            # Make graph - edges must be a list, and nodes must be a dic of
            # format {nodeindex: [x,y]}
            k = kn.KGraph(karstnet_edges, karstnet_nodes, verbose=False)
            metrics = k.characterize_graph(verbose)

            # Concatenate dataframes
            df_ = pd.DataFrame(metrics, index=[i])
            df_metrics = pd.concat([df_metrics, df_])

        return df_metrics

    @requires_karstnet()
    def compare_metrics(self,
                        dataframe: Union[DataFrame, Series, Styler],
                        ) -> Styler:
        """
        TODO

        Parameters
        ----------
        dataframe : Union[DataFrame, Series, Styler]
            Data to compare with karstnet metrics.

        Returns
        -------
        df_metrics : Styler
            
        Examples
        --------
        >>> app = pk.pykasso()
        >>> TODO
        >>> df_metrics = app.analyzer.compute_metrics()
        >>> app.analyzer.compare_metrics(df_metrics)
        """
        ### Convert pandas Series in DataFrame
        if isinstance(dataframe, Series):
            dataframe = dataframe.to_frame().T

        ### Define the text coloring function
        # Green if inside [V_min, V_max]
        # Orange if outside
        def _bg_color(x, min_val, max_val):
            if pd.isnull(x):
                return 'color: grey'
            elif (x < min_val) or (x > max_val):
                return 'color: #FF8C00'
            else:
                return 'color: #00FF00'

        # Iterate in the dataframe columns
        df_metrics = dataframe.style
        for column_name in dataframe:
            kwargs = {
                'min_val': self.stats[column_name]['min'],
                'max_val': self.stats[column_name]['max'],
                'subset': [column_name]
            }
            df_metrics = df_metrics.applymap(_bg_color, **kwargs)

        return df_metrics

    def compute_stats_on_networks(self,
                                  numpy_algorithm: str = 'mean',
                                  numpy_parameters: dict = {},
                                  ) -> np.ndarray:
        """
        Compute selected algorithm on the whole set of computed karst conduit
        networks.
        
        Parameters
        ----------
        numpy_algorithm : str, optional
            numpy algorithm, 'mean' by default.
            More details : https://numpy.org/doc/stable/reference/routines.statistics.html
        numpy_parameters : dict, optional
            Parameters of the selected algorithm, ``{}`` by default.
        
        Returns
        -------
        out : np.ndarray
        
        Examples
        --------
        >>> app = pk.pykasso()
        >>> TODO
        >>> df_metrics = app.analyzer.compute_metrics()
        >>> karst_std = app.analyzer.compute_stats_on_networks('std')
        """
        
        # For each simulation, retrieve data and store it
        karst_map = []
        for data in self.project:
            karst_map.append(data['maps']['karst'][-1].copy())
        
        # Retrieve algorithm
        try:
            numpy_func = getattr(np, numpy_algorithm)
        except ValueError:
            msg = "Asked algorithm is not valid."
            raise ValueError(msg)
        
        # Compute
        numpy_parameters.pop('axis', None)
        out = numpy_func(karst_map, axis=0, **numpy_parameters)
        
        return out
