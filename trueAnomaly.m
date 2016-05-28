function f = trueAnomaly(M,e)
% function f = trueAnomaly(M,e)
% Computes true anomaly for an elliptical orbit
% Required inputs:
%   M: mean anomaly [deg]
%   e: eccentricity [-]
% Outputs:
%   f: true anomaly [deg]

M = M*pi/180;

% Find eccentric anomaly
E = M;
E_err = 1;
E_tol = 1e-10;
while E_err > E_tol
    E_new = E + (M - E + e*sin(E))/(1 - e*cos(E));
    E_err = abs(E_new - E);
    E = E_new;
end

% Compute true anomaly
f = 2*atan(sqrt((1 + e)/(1 - e))*tan(E/2));
f = f*180/pi;
if f < 0
    f = f + 360;
end
