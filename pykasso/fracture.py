class Fracture():
    """
    Class modeling fractures as objects.
    """
    
    def __init__(self, ID, family, position, radius, orientation, dip, normal):
        """
        Create a fracture according to the parameters.
        
        Parameters
        ----------
        ID : int
            Fracture id.
        family : int
            Fracture family id.
        position : array
            Position of the center of the fracture [x, y, z].
        radius : float
            Radius of the fracture.
        orientation : float
            Orientation of the fracture.
        dip : float
            Dip of the fracture.
        normal : array
            Normal vector of the fracture [a, b, c].

        Examples
        --------
        >>> frac = pk.Fracture(0, 0, [10,10,10], 5, 180, 90, [1,0,0])
        """
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
        """
        Return the ID of the fracture.

        Examples
        --------
        >>> i = frac.get_ID()
        """
        return self.ID

    def get_family(self):
        """
        Return the family ID of the fracture.

        Examples
        --------
        >>> family = frac.get_family()
        """
        return self.family

    def get_position(self):
        """
        Return an array with the (x, y, z) coordinates of the center of the fracture.

        Examples
        --------
        >>> x, y, z = frac.get_position()
        """
        return self.position

    def get_radius(self):
        """
        Return the radius of the fracture.

        Examples
        --------
        >>> rad = frac.get_radius()
        """
        return self.radius

    def get_orientation(self):
        """
        Return the orientation of the fracture.

        Examples
        --------
        >>> orien = frac.get_orientation()
        """
        return self.orientation
        
    def get_dip(self):
        """
        Return the dip of the fracture.

        Examples
        --------
        >>> dip = frac.get_dip()
        """
        return self.dip

    def get_normal(self):
        """
        Return an array with the (a, b, c) vector component of the normal of the fracture.

        Examples
        --------
        >>> a, b, c = frac.get_normal()
        """
        return self.normal
