function [D,R,A] = sphericalTriangle(P1,P2,P3,outer,units)
% function [D,R,A] = sphericalTriangle(P1,P2,P3,outer,units)
% Returns the angular distances, rotation angles and total area of a
% spherical triangle defined by three points on the surface of a sphere.
% 
% Inputs:
%   - P1, P2, P3: the three points, each of which has to be a 1x2 vector,
%   with the first element indicating azimuth and the second one elevation.
%   - outer: a flag indicating whether the rotation angles and the area
%   should be computed for a large triangle (i.e. area > 2*pi steradian).
%   Can be either 1 (true) or 0 (false). The default value is 0.
%   - units: a string indicating the units in which P1, P2 and P3 are
%   provided. Can be either 'deg' or 'rad'. The default value is 'deg'.
% 
% Outputs:
%   - D: a 3x1 vector containing the angular distances between P1-P2,
%   P1-P3 and P2-P3, in the same units in which P1, P2 and P3 are provided.
%   - R: a 3x1 vector containing the rotation angles around P1, P2 and P3,
%   in the same units in which P1, P2 and P3 are provided.
%   - A: the total area of the spherical triangle, in steradians.
% 
% Author: Aleix Pinardell
% Version: 1.1.3
% Date: 25 February 2015
%

if nargin < 3
    error('You need to specify three points.');
end
if nargin < 4
    outer = 0;
end
if nargin < 5
    units = 'deg';
else
    units = validatestring(units,{'deg','rad'},'sphericalTriangle');
end
rightSize = [1 2];
if ~isequal(size(P1),rightSize) || ~isequal(size(P2),rightSize) || ...
        ~isequal(size(P3),rightSize)
    error('Wrong size for some input points. The right size is 1x2.');
end
if strcmp(units,'rad')
    P1 = rad2deg(P1);
    P2 = rad2deg(P2);
    P3 = rad2deg(P3);
end

P = [P1; P2; P3];
D = zeros(3,1);
R = zeros(3,1);
for i = 1:3
    for j = (i+1):3
        k = 6 - i - j;
        sinpart = sind(P(i,2))*sind(P(j,2));
        cospart = cosd(P(i,2))*cosd(P(j,2))*cosd(P(i,1) - P(j,1));
        D(k) = acosd(sinpart + cospart);
    end
end
for i = 1:3
    for j = (i+1):3
        k = 6 - i - j;
        cospart = cosd(D(k)) - cosd(D(j))*cosd(D(i));
        sinpart = sind(D(j))*sind(D(i));
        R(k) = acosd(cospart/sinpart);
    end
end

if outer
    R = 360 - R;
end
D = flip(D);
A = sum(deg2rad(R)) - pi;

if strcmp(units,'rad')
    D = deg2rad(D);
    R = deg2rad(R);
end
