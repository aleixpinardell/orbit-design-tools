function [x,y,z,v_x,v_y,v_z] = kep2car(a,e,i,raan,omega,M,mu)
% function [x,y,z,v_x,v_y,v_z] = kep2car(a,e,i,raan,omega,M,mu)
% Transforms Kepler coordinates of elliptic orbit to Cartesian coordinates
% Required inputs:
%   a: semi-major axis [km]
%   e: eccentricity [-]
%   i: inclination [deg]
%   raan: right ascension of the ascending node [deg]
%   omega: argument of pericenter [deg]
%   M: mean anomaly [deg]
% Optional inputs:
%   mu: gravitational parameter [km^3*s^-2]. If not specified, Earth's
%       gravitational parameter (398600.441 km^3*s^-2) is used.
% Outputs:
%   x,y,z: position vector in Cartesian coordinates [km]
%   v_x,v_y,v_z: velocity vector in Cartesian coordinates [km/s]

if nargin < 6 % if less than 6 inputs, break
    error('At least 6 input arguments required.');
else
    if e >= 1 % if parabollic or hyperbolic orbit, break
        error('The eccentricity must be smaller than 1.');
    end
    
    if nargin == 6 % if mu is not specified, use Earth's
        constants; % load constants from file
        mu = mu_E;
    end
    
    % Get true anomaly
    f = trueAnomaly(M,e);

    % Convert angles in deg to rad
    i = i*pi/180;
    raan = raan*pi/180;
    omega = omega*pi/180;
    f = f*pi/180;
    
    % Compute needed parameters for the transformation
    r = a*(1 - e^2)/(1 + e*cos(f));
    rf = [r*cos(f); r*sin(f)];
    l1 = cos(raan)*cos(omega) - sin(raan)*sin(omega)*cos(i);
    l2 = -cos(raan)*sin(omega) - sin(raan)*cos(omega)*cos(i);
    m1 = sin(raan)*cos(omega) + cos(raan)*sin(omega)*cos(i);
    m2 = -sin(raan)*sin(omega) + cos(raan)*cos(omega)*cos(i);
    n1 = sin(omega)*sin(i);
    n2 = cos(omega)*sin(i);
    
    % Determine position vector
    X = [l1 l2; m1 m2; n1 n2]*rf;
    x = X(1);
    y = X(2);
    z = X(3);
    
    % Determine velocity vector
    H = sqrt(mu*a*(1 - e^2));
    v_x = mu/H*(-l1*sin(f) + l2*(e + cos(f)));
    v_y = mu/H*(-m1*sin(f) + m2*(e + cos(f)));
    v_z = mu/H*(-n1*sin(f) + n2*(e + cos(f)));    
end
