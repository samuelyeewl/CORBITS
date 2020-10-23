'''
Python interface for CORBITS library by Brakensiek+Ragozzine

corbits.py
'''
import ctypes
import itertools
import numpy as np
import os.path

# Constants, identical to CORBITS
PI = 3.14159265358979
RAD_TO_DEG = 180 / PI
DEG_TO_RAD = 1/ RAD_TO_DEG
SR_TO_AU = 0.0046491
DAYS_IN_YEAR = 365.2425

dll_name = "libcorbits.so"
abs_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dll_path = abs_dir + os.path.sep + dll_name
clib = ctypes.CDLL(dll_path)

def stellar_mass_from_radius(radius):
    '''
    Convert stellar radius to mass, based on a broken power law:
    $$
        M_\star &= R_\star^{0.8},\qquad  R_\star >= 1.0 R_\odot \\
                &= R_\star^{0.57},\qquad R_\star < 1.0 R_\odot.
    $$

    Args:
    -----
    radius : float
        Stellar radius in solar units

    Returns:
    --------
    mass : float
        Stellar mass in solar units
    '''
    return (radius ** 0.8) if (radius >= 1.0) else (radius ** 0.57)


class sci_value(ctypes.Structure):
    '''
    Structure to contain value and uncertainties.
    '''
    _fields_ = [("val", ctypes.c_double),
                ("pos_err", ctypes.c_double),
                ("neg_err", ctypes.c_double)]

    def __init__(self):
        self.val = 0.0
        self.pos_err = 0.0
        self.neg_err = 0.0

    def __repr__(self):
        return '{:.6e} +{:6e}/-{:.6e}'.format(
            self.val, self.pos_err, self.neg_err)


class input_orbit(ctypes.Structure):
    '''
    Input orbit structure.

    Internally, angles are stored in radians; times in days; distances in AU.

    Parameters:
    -----------
    a : semimajor axis, default 1.0 AU
    e : eccentricity,   default 0.0
    i : inclination,    default 0.0
    omega : argument of periastron, default 0.0
    Omega : longitude of ascending node, default 0.0
    P : period, default 1.0 year (cannot provide both P and a)
    r_star : stellar radius, default 1.0 solar radii
    r : planet radius, default 0.0 earth radii
    use : whether to use planet, default 1

    Arguments:
    ----------
    angle_unit : 'deg' or 'rad'
    time_unit : 'day' or 'year'
    '''
    _fields_ = [("a", ctypes.c_double),
                ("r_star", ctypes.c_double),
                ("r", ctypes.c_double),
                ("e", ctypes.c_double),
                ("P", ctypes.c_double),
                ("Omega", ctypes.c_double),
                ("omega", ctypes.c_double),
                ("i", ctypes.c_double),
                ("use", ctypes.c_int)]

    def __init__(self, a=1.0, e=0.0, i=0.0, omega=0.0, Omega=0.0,
                 P=None, r_star=1.0, r=0.0, use=1,
                 angle_unit='deg', time_unit='year'):
        m_star = stellar_mass_from_radius(r_star)

        if P is not None:
            # Convert period to semimajor axis
            P_yr = P if time_unit == 'year' else P / DAYS_IN_YEAR
            a = pow(P_yr * P_yr * m_star, 1/3)
        else:
            # Convert semimajor axis to  period
            P = pow(a**3 / m_star, 0.5) * DAYS_IN_YEAR

        # Convert angles
        if angle_unit == 'deg':
            i *= DEG_TO_RAD
            omega *= DEG_TO_RAD
            Omega *= DEG_TO_RAD

        # Convert stellar radius to AU
        r_star *= SR_TO_AU

        super(input_orbit, self).__init__(
            a, r_star, r, e, P, Omega, omega, i, use)

clib.prob_of_transits_input_orbit.restype = sci_value
clib.prob_of_transits_input_orbit.argtypes = [ctypes.c_int,
                                              ctypes.POINTER(input_orbit)]
clib.prob_of_transits_input_orbit_noerr.restype = ctypes.c_double
clib.prob_of_transits_input_orbit_noerr.argtypes = [ctypes.c_int,
                                                    ctypes.POINTER(input_orbit)]

# Wrapper for prob_of_transits_input_orbit
def transit_prob(orbits, return_error=True):
    '''
    Compute the transit probability given the provided orbits.

    Args:
    -----
    orbits : list of input_orbit structures
    return_error: bool
        Whether to return error bounds

    Returns:
    --------
    sci_value or float
    '''
    if not isinstance(orbits, list):
        orbits = [orbits]
    n_pl = len(orbits)
    input_arr = (input_orbit * n_pl)(*orbits)

    if return_error:
        return clib.prob_of_transits_input_orbit(n_pl, input_arr)
    return clib.prob_of_transits_input_orbit_noerr(n_pl, input_arr)



