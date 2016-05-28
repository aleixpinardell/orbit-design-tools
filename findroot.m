function root = findroot(x, y)
% function root = findroot(x, y)
% 
% Finds the root of y(x) = 0.
% 
% Input parameters:
%   - x: a vector containing the values of the independent variable.
%   - y: a vector containing the values of the dependent variable.
% 
% Output parameter:
%   - root: the root of y(x) = 0, obtained by linear interpolation between
%     the two elements of y that are closest to zero.
% 
% Author: Aleix Pinardell
% Version: 1.0.0
% Date: 20 March 2016
% 

[~, index] = min(abs(y));
if y(index) < 0
    i = index;
    j = index + 1;
else
    i = index - 1;
    j = index;
end
root = x(j) - (x(j) - x(i))*y(j)/(y(j) - y(i));
