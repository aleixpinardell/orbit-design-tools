function value = H(angle)
% function value = H(angle)
% Returns the hemispheric function of an (array of) angles in radians.
% Namely, if mod(angle,2*pi) < pi, it returns 1. Otherwise it returns -1.
% 
% Author: Aleix Pinardell
% Version: 1.0
% Date: 26 February 2016
%

value = (mod(angle,2*pi) < pi) - (mod(angle,2*pi) >= pi);
