# uncertain_linear_interpolation
Propagate uncertainties through linear interpolation over irregular 
grids.

## Notes:
- Linear piecewise barycentric interpolation in N>2 dimensions 
using Delauney triangulation
- Identical to `scipy.interpolate.LinearNDInterpolator` but 
incorporates propagation of uncertainties through the interpolation
- Much slower than `scipy.interpolate.LinearNDInterpolator` since 
`scipy` uses just-in-time compilation, while this code does not
- IMPORTANT: The uncertainty propagation assumes that the value at each 
point are completely independent. Think carefully about this before you
use it.


TODO: Still needs documentation finished.
