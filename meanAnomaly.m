function M = meanAnomaly(f,e)
% function M = meanAnomaly(f,e)
% Computes mean anomaly for an elliptical orbit
% Required inputs:
%   f: true anomaly [deg]
%   e: eccentricity [-]
% Outputs:
%   M: mean anomaly [deg]

f = f*pi/180;
E = 2*atan(sqrt((1 - e)/(1 + e))*tan(f/2)); % compute eccentric anomaly
M = E - e*sin(E); % compute mean anomaly
M = M*180/pi;
if M < 0
    M = M + 360;
end
