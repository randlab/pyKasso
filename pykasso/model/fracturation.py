"""
This module contains fonctions designed for fracture generation and
discretization.
"""

import numpy as np
import pandas as pd


def _calculate_normal(
    dip: np.ndarray,
    orientation: np.ndarray,
) -> list:
    """
    Calculate the normal vectors for a set of fractures based on their dip and
    orientation angles.
    """
    orientation = np.radians(orientation)
    dip = np.radians(dip)
    x = np.sin(dip) * np.cos(np.radians(90) - orientation)
    y = np.sin(dip) * np.sin(np.radians(90) - orientation)
    z = np.cos(dip)
    a = y
    b = -x
    c = -z
    normal = list(zip(a, b, c))
    return normal


####################
### VOXELIZATION ###
####################


def _float_eq(a, b, tolerance: float = 1e-5) -> np.ndarray:
    """
    Returns True if the difference between a and b is lower than tolerance.

    Parameters
    ----------
    a : _type_
        _description_
    b : _type_
        _description_
    tolerance : float, optional
        _description_, by default 1e-5

    Returns
    -------
    result : bool
    """
    if np.all(np.isnan(a)):
        if np.all(np.isnan(b)):
            return True
        else:
            return False
    return np.all(abs(a - b) < tolerance)


def _unit_intersect(n: np.ndarray, d: float) -> np.ndarray:
    """
    Computes the intersection between a unit circle and a line on a 2D plane.
    The unit circle is centered on the origin (x=0, y=0) and has a radius of 1.
    The line is defined by parameters 'n' and 'd'.

    Parameters
    ----------
    n : numpy.ndarray of size 2
        The normal vector perpendicular to the line.
        Warning : The norm of n must be 1.
    d : float
        Defines the position of the line.
        It is a signed distance to the origin x=0, y=0.
        It can be positive or negative.

    Returns
    -------
    result : numpy.ndarray with 4 floating point values [X1, X2, Y1, Y2]
        The coordinates of the two points of intersection when it exists.
        Returns 4 np.nan if there is no intersection.
    """
    nx, ny = n[0], n[1]  # For readibility

    if not _float_eq(np.linalg.norm(n), 1):
        print("WARNING - unitcintl function : norm of n must be equal to 1")

    if np.abs(d) > 1:  # Case with no intersection
        return np.array([np.nan, np.nan, np.nan, np.nan])

    sd = np.sqrt(1 - d**2)

    if ny != 0:
        X1 = d * nx + ny * sd
        X2 = d * nx - ny * sd
        Y1 = d * ny - nx * sd
        Y2 = d * ny + nx * sd
    else:
        X1 = X2 = d
        Y1 = sd
        Y2 = -sd

    return np.array([X1, X2, Y1, Y2])


def _disk_zplane_intersect(
    center: np.ndarray,
    n: np.ndarray,
    R: float,
    zi: float,
) -> np.ndarray:
    """
    Computes the intersection between a disk and a horizontal plane (constant
    z plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : numpy.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : numpy.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    zi : float
        Position of the horizontal plane along the z-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [x1, y1, x2, y2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (zi - zc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, y1, x2, y2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[2] * (zc - zi) / tau / np.sqrt(1 - n[2]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, y1, x2, y2 = np.nan, np.nan, np.nan, np.nan

        else:
            n2 = n.copy()[0:2]  # Projection on horizontal plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            y1 = tau * xint[2] + yc
            y2 = tau * xint[3] + yc

    return np.array([x1, y1, x2, y2])


def _disk_xplane_intersect(
    center: np.ndarray,
    n: np.ndarray,
    R: float,
    xi: float,
) -> np.ndarray:
    """
    Computes the intersection between a disk and a vertical plane (constant x
    plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : np.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : np.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    xi : float
        Position of the vertical plane along the x-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [y1, z1, y2, z2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (xi - xc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        y1, z1, y2, z2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[0] * (xc - xi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            y1, z1, y2, z2 = np.nan, np.nan, np.nan, np.nan

        else:
            n2 = np.array([n[1], n[2]])  # Projection on vertical x plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            y1 = tau * xint[0] + yc
            y2 = tau * xint[1] + yc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([y1, z1, y2, z2])


def _disk_yplane_intersect(
    center: np.ndarray,
    n: np.ndarray,
    R: float,
    yi: float,
) -> np.ndarray:
    """
    Computes the intersection between a disk and a vertical plane (constant y
    plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : np.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : np.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    yi : float
        Position of the vertical plane along the y-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [x1, z1, x2, z2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (yi - yc)**2
    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, z1, x2, z2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[0] * (yc - yi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, z1, x2, z2 = np.nan, np.nan, np.nan, np.nan

        else:
            n2 = np.array([n[0], n[2]])  # Projection on vertical x plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([x1, z1, x2, z2])


def _rst2d(m: np.ndarray, xs: int, xe: int, ys: int, ye: int):
    """
    Rasterizes a line on a 2D plane.

    Parameters
    ----------
    m : numpy.ndarray of dim 2
    xs, xe, ys, ye : int
        Indices of positions on the grid m.
        Starting and ending location of the line.

    Returns
    -------
    m : numpy.ndarray of dim 2
        With value 1 where the grid is touched by the line.
    """

    nx, ny = m.shape

    if xe < xs:  # Ensuires dx always positive
        xe, xs = xs, xe
        ye, ys = ys, ye

    dx = xe - xs
    dy = ye - ys

    if xs >= nx:  # The line is entirely out of the grid
        return

    if (dx + np.abs(dy)) == 0:  # Case with a single pixel
        if (ys >= 0) and (ys < ny):
            m[xs, ys] = 1
        return

    if dx >= abs(dy):  # line with greater horizontal than vertical extension
        i = np.arange(max(0, xs), min(xe, nx))
        j = np.int32(np.round((i - xs) * dy / dx + ys))

        indx_ok = (j >= 0) & (j < ny)  # To crop the part that is out of grid
        i = i[indx_ok]
        j = j[indx_ok]

    else:
        if ye < ys:  # to ensure that arange works properly
            xe, xs = xs, xe
            ye, ys = ys, ye

        if ys > ny:
            return

        j = np.arange(max(ys, 0), min(ye, ny))
        i = np.int32(np.round((j - ys) * dx / dy + xs))

        indx_ok = (i >= 0) & (i < nx)  # To crop the part that is out of grid
        i = i[indx_ok]
        j = j[indx_ok]

    m[i, j] = 1
    return


def voxelize_fractures(
    grid,
    fractures: pd.DataFrame,
) -> np.ndarray:
    """
    Rasterizes a set of fractures on a 3D grid.
    
    Parameters
    ----------
    grid : _type_
        _description_
    fractures : _type_
        _description_

    Returns
    -------
    raster_fractures : numpy.ndarray of dim 3
        With value 1 where the grid is touched by a fracture.
    """

    dx, dy, dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    x0, y0, z0 = grid.x0, grid.y0, grid.z0

    # Creates empty array to store raster of fracture indicators
    raster_fractures = np.zeros((nx, ny, nz), dtype=np.int_)

    # Loop over the fractures
    for (i, f) in fractures.iterrows():

        # Get fracture geometry
        xc, yc, zc = f['x'], f['y'], f['z']
        nfx, nfy, nfz = f['normal']
        nfx, nfy, nfz = float(nfx), float(nfy), float(nfz)
        n = [nfx, nfy, nfz]
        R = f['radius']

        # Test orientation of the fractures
        if nfz**2 <= (nfx**2 + nfy**2):  # subvertical case
            # Vertical index of the center of the fracture
            kc = ((zc - z0) / dz)
            kc = int(kc)

            # Projected vertical extension
            vz = R * np.sqrt(1 - n[2]**2)

            # Vertical extension in number of cells
            dk = np.floor(vz / dz).astype(int)

            # Loop over the indices of horizontal levels
            for k in range(max(kc - dk, 0), min(kc + dk + 2, nz)):

                # Corresponding z value
                zi = z0 + k * dz + dz / 2

                # Computes the points of intersection of the fracture with
                # horizontal plane
                intersect = _disk_zplane_intersect([xc, yc, zc], n, R, zi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from x and y coordinates
                    i1, j1 = grid.get_i(intersect[0]), grid.get_j(intersect[1])
                    i2, j2 = grid.get_i(intersect[2]), grid.get_j(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:, :, k], i1, i2, j1, j2)

        elif nfx**2 < (nfz**2 + nfy**2):  # subhorizontal case 1
            # Horizontal x index of the center of the fracture
            ic = ((xc - x0) / dx)
            ic = int(ic)

            # Projected x extension
            vx = R * np.sqrt(1 - n[0]**2)
            # x extension in number of cells
            di = np.floor(vx / dx).astype(int)

            # Loop over the indices of horizontal levels
            for i in range(max(ic - di, 0), min(ic + di + 2, nx)):

                # Corresponding x value
                xi = x0 + i * dx + dx / 2

                # Computes the points of intersection of the fracture with
                # vertical x plane
                intersect = _disk_xplane_intersect([xc, yc, zc], n, R, xi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from y and z coordinates
                    j1, k1 = grid.get_j(intersect[0]), grid.get_k(intersect[1])
                    j2, k2 = grid.get_j(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[i, :, :], j1, j2, k1, k2)

        else:  # This may not be necessary
            # Horizontal y index of the center of the fracture
            jc = ((yc - y0) / dx)
            jc = int(jc)

            # Projected y extension
            vy = R * np.sqrt(1 - n[1]**2)
            # x extension in number of cells
            dj = np.floor(vy / dy).astype(int)

            # Loop over the indices of horizontal levels
            for j in range(max(jc - dj, 0), min(jc + dj + 2, nx)):

                # Corresponding x value
                yi = y0 + j * dy + dy / 2

                # Computes the points of intersection of the fracture with
                # vertical x plane
                intersect = _disk_yplane_intersect([xc, yc, zc], n, R, yi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from y and z coordinates
                    i1, k1 = grid.get_i(intersect[0]), grid.get_k(intersect[1])
                    i2, k2 = grid.get_i(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:, j, :], i1, i2, k1, k2)

        # raster_frac = np.swapaxes(raster_fractures,0,2)
    return raster_fractures
