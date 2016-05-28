function [TOFs, gamma1s] = alltofs(r1, r2, k2, Psi, N, stepSize, printProgress)
% function [TOFs, gamma1s] = alltofs(r1, r2, k2, Psi, N, stepSize, printProgress)
%
% Returns a vector of time of flights (TOFs) for an interplanetary flight
% modelled with an exposin, and the range of allowed initial flight path
% angles (gamma1s).
%
% Input parameters:
%   - r1: departure distance to central body [AU]
%   - r2: arrival distance to central body [AU]
%   - k2: winding parameter [-]
%   - Psi: true anomaly difference (<2*pi) [rad]
%   - N: number of complete revolutions [-]
%   - stepSize (optional): the step-size between two consecutive flight
%     path angles [rad]. The default is 0.01.
%   - printProgress (optional): a flag indicating whether the progress 
%     should be printed in the command window. The default is false.
% 
% Output parameters:
%   - TOFs: a vector containing the total time of flights [s].
%   - gamma1s: a vector containing the range of allowed initial flight path
%     angles [rad].
% 
% Requires:
%   - tof.m (the call to tof is done with no specified steps parameter)
% 
% See also tof
% 
% Author: Aleix Pinardell
% Version: 1.1
% Date: 20 March 2016
% 

if nargin < 7
    printProgress = 0;
end
if nargin < 6 || isempty(stepSize)
    stepSize = 0.01; % rad
end

theta_ = Psi + 2*pi*N;

% Compute bounds for gamma1
gamma1_term = -log(r1/r2)*cot(k2*theta_/2);
gamma1_Delta = 2*(1 - cos(k2*theta_))/k2^4 - log(r1/r2)^2;
gamma1_min = atan(k2/2*(gamma1_term - sqrt(gamma1_Delta)));
gamma1_max = atan(k2/2*(gamma1_term + sqrt(gamma1_Delta)));

% 
steps = ceil((gamma1_max - gamma1_min)/stepSize);
gamma1s = linspace(gamma1_min, gamma1_max, steps);
TOFs = zeros(size(gamma1s));
percent = -100;
for i = 1:length(gamma1s)
    TOFs(i) = tof(r1, r2, k2, Psi, N, gamma1s(i));
    if printProgress
        newPercent = floor(i/length(gamma1s)*100);
        if newPercent - percent >= 10
            fprintf('%g%%',newPercent);
            if newPercent < 100
                fprintf(' ... ');
            end
            percent = newPercent;
        end
    end
end
if printProgress
    fprintf('\n');
end

