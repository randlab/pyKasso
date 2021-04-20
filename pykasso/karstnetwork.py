class KarstNetwork():
    """
    A class for storing a calculated karst network.
	
	Parameters
    ----------
    maps : int
        -
    points : int
        -
    network : list
        -
    stats : float
        -
    settings : float
        -
    """
    def __init__(self, maps, points, network, stats, settings):
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings