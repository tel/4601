
# Fake data

from scipy.optimize import fmin
import numpy as np

def polarResolution(dx, w, h):

    def ellipseRadius(th):
        """
        For an ellipse of width w and height h, compute the radius at
        theta = th radians counter-clockwise away from a point on the
        major axis.

        In other words, rotate the ellipse so that it lies flat and
        then find the equation of the radius in polar form.
        """
        a = max([w, h])
        b = min([w, h])
        denom = sqrt((a * sin(th))**2 + (b * cos(th))**2)
        numer = a * b
        return numer/denom

    def maxError(dth):
        """
        Given a displacement in theta of dth, find the maximal error
        from the true radius.

        Maximal error will occur at the 'pointy' end of the ellipse,
        so we just subtract from the radius at theta = 0.
        """
        return 2*abs(ellipseRadius(0) - ellipseRadius(dth))

    # Find a maximizing estimate of what dth should be to get the
    # desired rectangular error.
    return mod(fmin(lambda dth: abs(dx - maxError(dth)), 0,
                    disp = False),
               pi)[0]

# sphere :: radius -> bounds -> ((z, th) -> distance)
def mkSphere(radius = 3, bounds = 6):
    """
    Given the radius of a sphere and the bounds of the measurement
    beam, return a function which computes the distance from one edge
    to the shadow of a sphere for some given theta and height.
    """
    def sphere(z, th):
        meas = bounds/2 - real(sqrt(2*radius*z - z**2))
        if isreal(meas):
            return meas
        else:
            return bounds
    return sphere

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

# scan :: (z, z) -> (dx, dz) -> ((z, th) -> distance) -> [(z, th, distance)]
def scan(zrange, rez, fn, w = 4, h = 4*1.15):

    # Compute ranges and resolutions
    zmin, zmax = zrange
    dx, dz = rez
    dth = polarResolution(dx, w, h)

    # Create sample positions
    zs  = linspace(zmin, zmax, num = 1/dz)
    ths = linspace(0, 2*pi, 1/dth)
    samps = cartesian((zs, ths))

    return  apply_along_axis(lambda x: r_[x, fn(*x)], 1, samps)

fakesphere = lambda: scan([0, 6], [0.002, 0.002], mkSphere(3, 6))
tinysphere = lambda: scan([0, 1], [0.2, 0.2], mkSphere(1, 2))

# reconstructP :: bounds -> [(z, th, distance)] -> [(z, th, r)]
def reconstructP(bounds, samples):
    def renorm(dims):
        z, th, dist = dims
        return array([z, th, bounds/2-dist])
    return apply_along_axis(renorm, 1, samples)

# reconstructR :: [(z, th, r)] -> [(z, x, y)]
def reconstructR(samples):
    def renorm(dims):
        z, th, r = dims
        return array([z, r*sin(th), r*cos(th)])
    return apply_along_axis(renorm, 1, samples)
