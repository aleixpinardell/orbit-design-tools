function vout = rowVectorOfSize(vin, n)
% function vout = rowVectorOfSize(vin, n)
% - If vin is an empty vector, returns a vector of size 1xn with all its elements set to zero.
% - If vin is a scalar, returns a vector of size 1xn with all its elements set to vin.
% - If vin is a vector of size 1xn, returns vin.
% - If vin is a vector of size nx1, returns the transpose of vin.
% - Otherwise, the function throws an error.
% Author: Aleix Pinardell
% Version: 1.0
% Date: 26 November 2015

s = size(vin);
o = ones(1,n);
vout = vin;
if isempty(vin)
    vout = 0*o;
elseif s(1) ~= 1 || s(2) ~= n
    if length(vin) == 1
        vout = vin*o;
    elseif s(1) == n && s(2) == 1
        vout = vin';
    else
        error('Inconsistent sizes.');
    end
end
