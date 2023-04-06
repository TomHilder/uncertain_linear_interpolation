import numpy as np
import scipy.spatial.distance as dist

from scipy.spatial import Delaunay


def find_simplex(points, targets):
    """Find the simplex of each target point that contains it from the points array using Delaunay triangulation.

    Args:
        points (numpy.ndarray): A numpy array representing a set of N-dimensional points with shape n x N.
        targets (numpy.ndarray): A numpy array representing another set of N-dimensional points with shape m x N.

    Returns:
        numpy.ndarray: A numpy array representing the indices of the vertices of the simplex for each target point. The output has shape m x (N+1), where each row is the index of the vertices of the simplex for that target point.
    """

    # Compute the Delaunay triangulation of the points
    tri = Delaunay(points)

    # Find the indices of the simplices that contain each target point
    simplices = tri.find_simplex(targets)

    # Return the indices of the vertices of the simplices for each target point
    return tri.simplices[simplices], simplices


def barycentric_coordinates(point, simplex):
    """Find the barycentric coordinates of a point in a simplex in N dimensions.

    Args:
        point (numpy.ndarray): A numpy array representing the point whose barycentric coordinates are to be found.
        simplex (numpy.ndarray): A numpy array representing the simplex in which the point lies. The simplex is assumed to be an (N+1) x N numpy array, where N is the number of dimensions.

    Returns:
        numpy.ndarray: A numpy array representing the barycentric coordinates of the point in the simplex.
    """

    # Compute the matrix A = [v1 - vN, v2 - vN, ..., vN-1 - vN]
    A = simplex[:-1] - simplex[-1]

    # Compute the vector b = p - vN
    b = point - simplex[-1]

    # Compute the barycentric coordinates by solving the system Ax = b
    x = np.linalg.solve(A.T, b)

    # Add a final coordinate of 1-x1-x2-...-xN to the end of x to obtain the full set of barycentric coordinates
    return np.append(x, 1 - np.sum(x))


def cartesian_coordinates(barycentric, simplex):
    """Find the Cartesian coordinates of a point given its barycentric coordinates and the simplex in which it lies.

    Args:
        barycentric (numpy.ndarray): A numpy array representing the barycentric coordinates of the point.
        simplex (numpy.ndarray): A numpy array representing the simplex in which the point lies. The simplex is assumed to be an (N+1) x N numpy array, where N is the number of dimensions.

    Returns:
        numpy.ndarray: A numpy array representing the Cartesian coordinates of the point.
    """

    # Multiply the barycentric coordinates by the simplex vertices, and sum the result along the rows to obtain the Cartesian coordinates
    return np.sum(barycentric.reshape(-1, 1) * simplex, axis=0)



def find_all_barycentric_coords(points, targets):
    
    # find simplices
    simplices, inds = find_simplex(points, targets)
    
    # loop over all targets
    bary_coords = []
    for i in range(len(inds)):
        
        # get barycentric coords for point in simplex
        bary = barycentric_coordinates(targets[i,:], points[simplices[i],:])


        if inds[i] > 0:
            bary_coords.append(bary)
        else:
            nans_array    = np.empty(bary.shape)
            nans_array[:] = np.nan
            bary_coords.append(nans_array)
            
    return np.array(bary_coords), simplices