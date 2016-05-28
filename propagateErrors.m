function [rss, errors] = propagateErrors(alt,el,az,lat,sourceErrors,propMode)
% [rss, errors] = propagateErrors(alt,el,az,lat,sourceErrors,propMode)
% 
% Inputs:
%   * alt: altitude of the spacecraft [km]
%   * el: elevation angle of the spacecraft as seen from the ground [deg]
%   * az: azimuth angle of the spacecraft as seen from the ground [deg]
%   * lat: latitude of the target [deg]
%   * sourceErrors: 1x7 vector containing the following source errors:
%       - azimuth [deg]
%       - nadir angle [deg]
%       - in-track position [km]
%       - cross-track position [km]
%       - radial position [km]
%       - altitude [km]
%       - spacecraft clock [s]
%   * propMode: either 'mapping' or 'pointing' (default is 'mapping')
% Outputs:
%   * rss: root sum square of the propagated error [in km for mapping
%     errors, in deg for pointing errors]
%   * errors: 1x7 array containing the individual propagated errors
% 
% Author: Aleix Pinardell
% Version: 1.0.1
% Date: 17 February 2015

if nargin < 6
    propMode = 'mapping';
end
sourceErrors(1:2) = pi/180*sourceErrors(1:2);

% Load Earth-related parameters
constants;

% Geometry
el = deg2rad(el);
az = deg2rad(az);
lat = deg2rad(lat);
R_T = R_E; % distance from Earth's center to the target [km]
R_S = R_E + alt; % distance from Earth's center to the satellite [km]
rho = asin(R_E/R_S); % angular radius of the Earth [rad]
V_e = 2*pi*R_E/(3600*24); % Earth rotation velocity at the equator [km/s]
nadir = asin(cos(el)*sin(rho)); % [rad]
central = pi/2 - el - nadir; % [rad]
D = R_E*sin(central)/sin(nadir); % target distance [km]
H = asin(sin(central)*sin(az));
G = asin(sin(central)*cos(az));
Y_I = acos(cos(az)*sin(nadir));
Y_C = acos(sin(az)*sin(nadir));
J = acos(cos(az-pi/2)*cos(el));


% Error propagation
if strcmp(propMode,'mapping')
    propagationVector = [D*sin(nadir);
                         D/sin(el);
                         R_T/R_S*cos(H);
                         R_T/R_S*cos(G);
                         sin(nadir)/sin(el);
                         1/tan(el);
                         V_e*cos(lat)];
elseif strcmp(propMode,'pointing')
    propagationVector = [sin(nadir);
                         1;
                         sin(Y_I)/D;
                         sin(Y_C)/D;
                         sin(nadir)/D;
                         0;
                         V_e/D*cos(lat)*sin(J)];
    propagationVector = propagationVector*180/pi;
else
    error('The value provided for ''propMode'' is not valid.');
end
errors = sourceErrors.*propagationVector';
rss = sqrt(errors*errors');
