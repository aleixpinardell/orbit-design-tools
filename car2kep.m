function [a,e,i,raan,omega,M] = car2kep(X,V,mu)
% function [a,e,i,raan,omega,M] = car2kep(X,V,mu)
% Transforms Cartesian coordinates of elliptic orbit to Kepler coordinates
% Inputs:
%   X = [x,y,z,v_x,v_y,v_z] state vector in Cartesian coordinates [km...km/s]
%   If V = [v_x,v_y,v_z] is specified, then X = [x,y,z]
% Optional inputs:
%   mu: gravitational parameter [km^3*s^-2]. If not specified, Earth's
%       gravitational parameter (398600.441 km^3*s^-2) is used.
% Outputs:
%   a: semi-major axis [km]
%   e: eccentricity [-]
%   i: inclination [deg]
%   raan: right ascension of the ascending node [deg]
%   omega: argument of pericenter [deg]
%   M: mean anomaly [deg]

if nargin == 0 % if no inputs, break
    error('At least 1 input argument required.');
else    
    if nargin < 3 % if mu is not specified, use Earth's
        constants; % load constants from file
        mu = mu_E;
    end
    if nargin == 1 % if X is a state-vector and V is not specified
        if length(X) == 6
            V = X(4:6);
            X = X(1:3);
        else
            error('The state-vector has to have 6 components.');
        end
    else % check X and V have 3 components each
        if length(X) ~= 3 || length(V) ~= 3
            error('The position and velocity vectors have to have 3 components each.');
        end
    end
    
    r = norm(X);
    v = norm(V);
    h = cross(X,V);
    N = cross([0;0;1],h);
    
    a = 1/(2/r - v^2/mu);
    e = 1/mu*cross(V,h) - X/r;
    i = acos(h(3)/norm(h));
    Nxy = sqrt(N(1)^2 + N(2)^2);
    raan = atan2(N(2)/Nxy,N(1)/Nxy);
    
    NN = N/norm(N);
    ee = e/norm(e);
    omega = sign(dot(cross(NN,e),h))*acos(dot(ee,NN));
    f = sign(dot(cross(e,X),h))*acos(dot(X/r,ee));
    
    e = norm(e);
    i = i*180/pi;
    if isnan(raan) % If the raan is not defined, assign (arbitrary) 0 value
        raan = 0;
    else
        raan = raan*180/pi;
        if raan < 0
            raan = raan + 360;
        end
    end
    if isnan(omega) % If omega is not defined, assign (arbitrary) 0 value
        omega = 0;
    else
        omega = omega*180/pi;
        if omega < 0
            omega = omega + 360;
        end
    end
    f = f*180/pi;
    if f < 0
        f = f + 360;
    end
    M = meanAnomaly(f,e);
end
