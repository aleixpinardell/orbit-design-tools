function [A,c,B] = sideAngleSide(a,C,b,second,rad)
% function function [A,c,B] = sideAngleSide(a,C,b,second,rad)
% Solves a spherical triangle defined by side-angle-side. All the rotation
% angles are defined by the right-hand rule, inside-out.
% 
% Inputs:
%   - a: the angular distance of the first side.
%   - C: the rotation angle from the first side to the second side.
%   - b: the angular distance of the second side.
%   - second (optional): a flag specifying whether the solution to be
%     returned is the first one (0) or the second one (1). Default is 0.
%   - rad (optional): a flag specifying whether the angles are provided and
%     returned in radians (1) or degrees (0). The default value is 0.
% 
% Outputs:
%   - A: the rotation angle from the second side to the third side.
%   - c: the angular distance of the third side.
%   - B: the rotation angle from the thrid side to the first side.
% 
% Author: Aleix Pinardell
% Version: 1.0.2
% Date: 29 February 2016
%

if nargin < 5
    rad = 0;
    if nargin < 4
        second = 0;
    end
end
if rad
    a = rad2deg(a);
    C = rad2deg(C);
    b = rad2deg(b);
end

c = acos2d(validate(cosd(a)*cosd(b) + sind(a)*sind(b)*cosd(C),-1,1), Hd(C));
A = acos2d(validate((cosd(a)-cosd(b)*cosd(c))/(sind(b)*sind(c)),-1,1), Hd(a));
B = acos2d(validate((cosd(b)-cosd(a)*cosd(c))/(sind(a)*sind(c)),-1,1), Hd(b));

if second
    c = 360 - c;
    A = mod(A + 180, 360);
    B = mod(B + 180, 360);
end
if rad
    A = deg2rad(A);
    c = deg2rad(c);
    B = deg2rad(B);
end
