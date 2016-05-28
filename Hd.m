function value = Hd(angle)
% function value = Hd(angle)
% Returns the hemispheric function of an (array of) angles in degrees.
% Namely, if mod(angle,360) < 180, it returns 1. Otherwise it returns -1.
% 
% Author: Aleix Pinardell
% Version: 1.0
% Date: 26 February 2016
%

value = (mod(angle,360) < 180) - (mod(angle,360) >= 180);
