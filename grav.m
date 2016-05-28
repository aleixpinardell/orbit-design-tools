function dXdt = grav(t,X)
% function dXdt = grav(t,X)
% Inputs:
%   t: time, irrelevant for this case.
%   X: (2n)x1 state vector, with first half of components indicating
%      position and last half of components indicating velocities.
% Output:
%   dXdt: time derivative of the state vector evaluated at X and t. Same size
%         as X. Only Earth's main gravitational acceleration is considered.
%
%  Author: Aleix Pinardell
%  Version: 1.2
%  Date: 4 December 2015

constants;
dXdt = zeros(size(X));
n = ceil(length(X)/2); % this is generally (but not necessarily) 3
dXdt(1:n) = X(n+1:2*n); % store velocity in first n components
r = X(1:n);
dXdt(n+1:2*n) = -mu_E/norm(r)^3*r; % store acceleration in last n components
