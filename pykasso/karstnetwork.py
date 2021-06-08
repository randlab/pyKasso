class KarstNetwork():
    """
    Class modeling the computed karst network and containing all the data used to get it.
    """

    def __init__(self, maps, points, network, stats, settings):
        """
        Creates an object and stores within all the results computed by SKS instance.

        Parameters
        ----------
        maps : array
            Data.
        points : array
            Inlets and oulets used for the simulation.
        network : array

        stats : array
            Statistics calculated for the karst network.
        settings : dict
            Dictionnary of the settings used for the simulation.
        """
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings
