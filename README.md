# Orbit design tools
MATLAB Tools for Mission Geometry &amp; Orbit Design

## General

### Numerical integration
#### `nint`
Integrates a state vector using the a numerical integration method.
#### `grav`
Simplest gravitational model to be used for integration.

### Optimisation
#### `newrap`
Finds the minimum of a function F(x) using the Newton-Raphson method.
#### `geneticAlgorithm`
Finds the optimum of a function f(x) using a genetic algorithm.

### Spherical trigonometry
#### `sphericalTriangle`
Returns the angular distances, rotation angles and total area of a spherical triangle defined by three points on the surface of a sphere.

### Other
#### `findroot`
Finds the root of `y(x) = 0`.

#### `rowVectorOfSize`
Returns an array of size 1xn

#### `validate`
Checks that a value is within the range `[lowerLimit - tolerance, upperLimit + tolerance]`. If it is, the function returns the input value constrained to the range `[lowerLimit, upperLimit]`. It it is not, the function throws an error.

## Astrodynamics

#### `constants.m`
Loads astrodynamics constants in the Workspace

#### `meanAnomaly`
Computes mean anomaly for an elliptical orbit
#### `trueAnomaly`
Computes true anomaly for an elliptical orbit

### Coordinate systems transformations
#### `car2kep`
Transforms Cartesian coordinates of elliptic orbit to Kepler coordinates
#### `kep2car`
Transforms Kepler coordinates of elliptic orbit to Cartesian coordinates

### Earth-repeat orbits
#### `earthRepeatOrbit`
Considers a satellite in an Earth repeat orbit completing `j` orbits every
`k` days. The numerical value of two of the following orbital parameters
must be specified:
- semi-major axis
- eccentricity
- inclination

The non-specified parameter must be introduced as an empty vector `[]`.
The function solves for the value of the non-specified parameter that
satisfies the requirements on the other two parameters, `j` and `k`.


### Satellite mapping and pointing errors
#### `propagateErrors`
Returns the propagated mapping and pointing errors for a given satellite


### Full-sky geomtry
#### `acos2`
Implements the function defined in Orbit & Constellation Design & Management, Appendix A.7.7. Output in radians.
#### `acos2d`
Implements the function defined in Orbit & Constellation Design & Management, Appendix A.7.7. Output in degrees.

#### `H`
Returns the hemispheric function of an (array of) angles in radians.
#### `Hd`
Returns the hemispheric function of an (array of) angles in degrees.

#### `sideAngleSide`
Solves a spherical triangle defined by side-angle-side. All the rotation angles are defined by the right-hand rule, inside-out.

#### `dualAxis`
Solves a dual-axis spiral problem.


### Interplanetary trajectories
#### `tof`
Returns the time of flight (TOF) of an interplanetary flight modelled with an exposin.

#### `alltofs`
Returns a vector of time of flights (TOFs) for an interplanetary flight modelled with an exposin, and the range of allowed initial flight path angles.
