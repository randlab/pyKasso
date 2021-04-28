class Fracture():
    """
    A class for modeling fractures as objects.

    Parameters
    ----------
    ID : int
        Fracture id.
    family : int
        Fracture family id.
    position : list
        Position of the center of the fracture [x, y, z].
    radius : float
        Radius of the fracture.
    orientation : float
        Orientation of the fracture.
    dip : float
        Dip of the fracture.
    normal : list
        Normal vector of the fracture [a, b, c].
    """

    def __init__(self, ID, family, position, radius, orientation, dip, normal):
        self.ID          = ID
        self.family      = family
        self.position    = position
        self.radius      = radius
        self.orientation = orientation
        self.dip         = dip
        self.normal      = normal

    def __repr__(self):
        return '[id:{}, fam.:{}, x:{}, y:{}, z:{}, rad:{}, or.:{}, dip.:{}, n:({},{},{})] \n'.format(self.ID,self.family, round(self.position[0],2), round(self.position[1],2), round(self.position[2],2), round(self.radius,2), round(self.orientation,2), round(self.dip,2), round(self.normal[0],2), round(self.normal[1],2), round(self.normal[2],2))

    def get_ID(self):
        return self.ID

    def get_family(self):
        return self.family

    def get_position(self):
        return self.position

    def get_radius(self):
        return self.radius

    def get_orientation(self):
        return self.orientation
        
    def get_dip(self):
        return self.dip

    def get_normal(self):
        return self.normal