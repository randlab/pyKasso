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
from pykasso._typing import Project, Series, DataFrame, Styler


def requires_karstnet():
    """
    If 'karstnet' package is not installed, return `ImportError` exception
    when a method requiring 'karstnet' is called.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            if not _has_karstnet:
                msg = ("karstnet is required to do this.")
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
    def __init__(self, project: Project) -> None:
        """
        Initialize the class.
        """
        # Load reference metrics for statistical karstic network analysis
        self.project = project
        self.k = None  # TODO : to remove ?
        self._load_statistics()
        
    def _load_statistics(self) -> None:
        """
        Set the reference metrics for statistical karstic network analysis
        More details here : https://github.com/karstnet/karstnet
        """
        package_location = self.project._pckg_paths['package_location']
        statistics_file_path = "/../_misc/statistics.xlsx"
        statistics_file_location = package_location + statistics_file_path
        self.statistics = pd.read_excel(statistics_file_location).describe()
        return None
    
    @requires_karstnet()
    def compute_metrics(self) -> DataFrame:
        """
        Computes the statistical metrics for each simulated karst network
        using the karstnet package.

        Returns
        -------
        df_metrics : pd.DataFrame
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
        >>> import pykasso.analysis as pka
        >>> analysis = pka.Analysis('examples/betteraz/')
        >>> metrics = analysis.compute_karstnet_metrics()
        """
        df_metrics = pd.DataFrame()

        # For each simulation, retrieves data and computes metrics
        for i, path in enumerate(self.project.simulations):
            # Retrieves data
            results = self.project._read_pickle(path + "results.pickle")
            karstnet_edges = list(results["vectors"]["edges"].values())
            karstnet_nodes = copy.deepcopy(results["vectors"]["nodes"])
            
            # Drops last item in list (the node type) for each dictionary entry
            for key, value in karstnet_nodes.items():
                value.pop()

            # Computes karstnet metrics
            # Makes graph - edges must be a list, and nodes must be a dic of
            # format {nodeindex: [x,y]}
            self.edges = karstnet_edges  # TODO - remove self
            self.nodes = karstnet_nodes  # TODO - remove self
            # TODO - remove 'self' ?
            self.k = kn.KGraph(karstnet_edges, karstnet_nodes)
            verbosity = 0
            metrics = self.k.characterize_graph(verbosity)

            # Concatenates dataframes
            df_ = pd.DataFrame(metrics, index=[i])
            df_metrics = pd.concat([df_metrics, df_])

        return df_metrics

    def compare_metrics(self, data: Union[Styler, DataFrame]) -> Styler:
        """_summary_

        Parameters
        ----------
        data : Series, DataFrame
            _description_

        Returns
        -------
        df_metrics :
            _description_
            
        Examples
        --------
        >>> import pykasso.analysis as pka
        >>> analysis = pka.Analysis('examples/betteraz/')
        >>> metrics = analysis.compute_karstnet_metrics()
        >>> analysis.compute_karstnet_metrics(metrics.mean())
        """
        # Converts pandas Series in DataFrame
        if isinstance(data, (pd.Series)):
            data = data.to_frame().T

        # Defines the text coloring function
        # Red if
        # Bright red if
        # Green if
        def _bg_color(x, min_val, max_val, mean_val, std_val):
            if (x < min_val) or (x > max_val):
                return 'color: red'
            else:
                inf = mean_val - std_val
                sup = mean_val + std_val
                if (x < inf) or (x > sup):
                    return 'color: #FFFF00'
                else:
                    return 'color: #00FF00'

        # Iterates in the dataframe columns
        df_metrics = data.style
        for column_name in data:
            kwargs = {
                'min_val': self.statistics[column_name]['min'],
                'max_val': self.statistics[column_name]['max'],
                'mean_val': self.statistics[column_name]['mean'],
                'std_val': self.statistics[column_name]['std'],
                'subset': [column_name]
            }
            df_metrics = df_metrics.applymap(_bg_color, **kwargs)

        return df_metrics

    def compute_stats_on_networks(self, algorithm: str = 'mean',
                                  np_options: dict = {}) -> np.ndarray:
        """
        TODO
        
        Computes the mean of all the simulations.
        
        Returns
        -------
        out : np.ndarray
            __doc__
        
        Examples
        --------
        >>> import pykasso.analysis as pka
        >>> analysis = pka.Analysis('examples/betteraz/')
        >>> karst_mean = analysis.compute_average_karstic_network()
        """
        # Initialization
        # nx = self.project.grid.nx
        # ny = self.project.grid.ny
        # nz = self.project.grid.nz
        karst_map = []
        # karst_map = np.zeros((nx, ny, nz))
        
        # For each simulation, retrieves data and sums it # TODO
        for path in self.project.simulations:
            results = self.project._read_pickle(path + 'results.pickle')
            # karst_map = np.add(karst_map, results['maps']['karst'][-1]).copy()
            karst_map.append(results['maps']['karst'][-1].copy())
        
        # Calculates
        try:
            numpy_func = getattr(np, algorithm)
        except ValueError:
            msg = "Asked algorithm is not valid."  # TODO
            raise ValueError(msg)
        
        np_options.pop('axis', None)
        out = numpy_func(karst_map, axis=0, **np_options)
        return out
