#include <math.h>
#include <stdlib.h>

/* miscelaneous mathematical and astronomical functions and constants */

/* --- CONSTANTS --- */

double const pi = 3.14159265358979;

double const rad_to_deg = 180 / pi;

double const eps = 2e-9; /* margin of error in floating-point calculations */

double const SRtoAU = .0046491; /* http://en.wikipedia.org/wiki/Solar_radius */

double const days_in_year = 365.2425;

/* wikipedia */
double const SR_to_km = 696342;

/* wikipedia */
double const JR_to_km = 69911;

double const solar_mass = 1.989e30; /* http://solarscience.msfc.nasa.gov/ */

double const SM_to_EM = 332946; /* http://en.wikipedia.org/wiki/Solar_Mass */

/* Murray and Correia: Keplerian Orbits and Dynamics of Exoplanets*/
double const G = 6.67260e-11; /* N (m/kg)^2 */

/* --- FUNCTION PROTOTYPES --- */

/* absolute value of the double x */
double abs_double(double x);

/* determines if two doubles are equil within an eps of error */
/* returns 1 if equal, 0 otherwise */
int equal(double a, double b);

/* returns the minimum of two doubles */
double min(double a, double b);

/* semimajor axis given mass of star and period of orbit */
/* mass is in solar masses; period is in days*/
/* From Kepler's Third Law, see <http://www.astro.lsa.umich.edu/undergrad/
   Labs/extrasolar_planets/pn_intro.html> */
double radius(double M_star, double P);

/* approximate mass of a star, given only its solar radius */
/* based on http://www2.astro.psu.edu/users/rbc/a534/lec18.pdf */
/* input: SR, output: SM */
double mass_from_radius (double R_star);

/* returns a random number from the uniform distribution for [0, range]*/
double rand_uniform (double range);

/* returns a random number from the Rayleigh distribution with width sigma */
double rand_Rayleigh (double sigma);