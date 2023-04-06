from .geometry import *

def linear_interpolation_ND(points, values, u_values, new_points):
    
    # get barycentric coordinates and simplices for all new points
    barycentric_coords, points_inds = find_all_barycentric_coords(points, new_points)
    
    # weighted sums over values at old points using barycentric coords gives values at new points 
    # (this is linear interpolation)
    new_values = np.sum(barycentric_coords * values[points_inds], axis=1)

    # now propagate uncertainties (yes, I wrote this whole code just to get this step)
    new_u_values = np.sqrt(
        np.sum(barycentric_coords**2 * u_values[points_inds]**2, axis=1)
    )
    
    return new_values, new_u_values