function value = validate(value, lowerLimit, upperLimit, tolerance)
% function value = validate(value, lowerLimit, upperLimit, tolerance)
% Checks that a value is within the range [lowerLimit - tolerance,
% upperLimit + tolerance]. If it is, the function returns the input value
% constrained to the range [lowerLimit, upperLimit]. It it is not, the
% function throws an error.
% 
% Inputs:
%   - value: the value to be checked.
%   - lowerLimit: the smallest acceptable value.
%   - upperLimit (optional): the largest acceptable value. If not
%     specified, it is assumed to be equal to lowerLimit.
%   - tolerance (optional): the value used to extend the allowable range.
%     Default value is 1e-5.
% Output:
%   - value: the value constrained to the range [lowerLimit, upperLimit].
% 
% Author: Aleix Pinardell
% Version: 1.0
% Date: 26 February 2016
%

if nargin < 4
    tolerance = 1e-5;
    if nargin < 3
        upperLimit = lowerLimit;
    end
end

if ~isnumeric([value lowerLimit upperLimit tolerance])
    error('The input arguments of this function must be numeric.');
end

for i = 1:length(value)
    min = lowerLimit - tolerance;
    max = upperLimit + tolerance;
    if value(i) < min || value(i) > max
        error('%g is not within the allowable range [%g, %g]',value(i),min,max);
    elseif value(i) < lowerLimit
        value(i) = lowerLimit;
    elseif value(i) > upperLimit
        value(i) = upperLimit;
    end
end
