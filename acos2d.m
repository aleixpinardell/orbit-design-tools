function angle = acos2d(value,h)
% function angle = acos2d(value)
% Implements the function defined in Orbit & Constellation Design &
% Management, Appendix A.7.7
% 
% Inputs:
%   - value: from -1 to 1.
%   - h: the value of the hemispheric function. Can be either 1 or -1.
% Output: angle (in deg).
% 
% Author: Aleix Pinardell
% Version: 1.0.1
% Date: 26 February 2016
%

angle = mod(h.*acosd(value),360);
