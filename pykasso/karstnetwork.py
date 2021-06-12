class KarstNetwork():
    """
    Class modeling the computed karst network and containing all the data used to get it.
    """

    def __init__(self, maps, points, network, stats, settings):
        """
        Creates an object and stores within all the results computed by the SKS instance.

        Parameters
        ----------
        maps : array
            Data.
        points : array
            Inlets and outlets used for the simulation.
        network : array
            Karst network.
        stats : array
            Statistics calculated for the karst network.
        settings : dict
            Dictionary of the settings used for the simulation.
        """
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings
