function [TOF, theta, s, c, tan_gamma, gamma, r, dt, t] = ...
    tof(r1, r2, k2, Psi, N, gamma1, steps)
% function [TOF, theta, s, c, tan_gamma, gamma, r, dt, t] = ...
%    tof(r1, r2, k2, Psi, N, gamma1, steps)
% 
% Returns the time of flight (TOF) of an interplanetary flight modelled
% with an exposin.
%
% Input parameters:
%   - r1: departure distance to central body [AU]
%   - r2: arrival distance to central body [AU]
%   - k2: winding parameter [-]
%   - Psi: true anomaly difference (<2*pi) [rad]
%   - N: number of complete revolutions [-]
%   - gamma1: initial flight path angle [rad]
%   - steps (optional): number of steps in which the integration is
%     divided. If not specified, the default value is initially set to 10.
%     Then, the program iterates while duplicating the number of steps
%     until the computed TOF differs less than 0.001% with respect to the
%     previous value.
%
% Ouput parameters:
%   - TOF: total time of flight [s]
%   - theta: vector with the values of the true anomaly at each step [rad]
%   - s: vector with the values of s at each step [rad]
%   - c: vector with the values of c at each step [rad]
%   - tan_gamma: vector with the values of the tangent of the flight path
%     angle at each step [rad]
%   - gamma: vector with the values of the flight path angle at each step
%     [rad]
%   - r: radial distance to central body at each step [km]
%   - dt: vector with the time increments at each step [s]
%   - t: vector with the accumulated time at each step [s]
% 
% Requires:
%   - constants.m
% 
% See also constants
% 
% Author: Aleix Pinardell
% Version: 1.0.2
% Date: 20 March 2016
% 

theta_ = Psi + 2*pi*N;

% Loads mu_S (gravitational parameter of the Sun) and AU (1 AU in km)
constants;

% Compute k1, phi and k0
k1_term = log(r1/r2) + tan(gamma1)/k2*sin(k2*theta_);
k1 = sign(k1_term)*sqrt((k1_term/(1 - cos(k2*theta_)))^2 + tan(gamma1)^2/k2^2);
phi = acos(tan(gamma1)/(k1*k2));
k0 = r1/(exp(k1*sin(phi)))*AU;

% Determine if the number of steps if fixed or should try converge TOF
fixedNumberOfSteps = 1;
if nargin < 7 || isempty(steps)
    fixedNumberOfSteps = 0;
    steps = 10;
end

TOF_prev = 0;
tol = 1e-5; % relative
change = 2*tol;
while change > tol
    % Integration steps
    theta = linspace(theta_/(2*steps), theta_*(2*steps-1)/(2*steps), steps);
    arg = k2*theta + phi;
    s = sin(arg);
    c = cos(arg);
    tan_gamma = k1*k2.*c;
    gamma = atan(tan_gamma);
    r = k0.*exp(k1.*s);

    % Find time of flight
    integrand = sqrt(r.^3.*(tan_gamma.^2 + k1.*k2^2.*s + 1)/mu_S);
    dt = integrand*theta_/steps;
    t = zeros(size(dt));
    for i = 1:length(t)
        t(i) = sum(dt(1:i));
    end
    TOF = t(end);
    if fixedNumberOfSteps
        break;
    else
        change = abs(TOF - TOF_prev)/TOF_prev;
        TOF_prev = TOF;
        steps = 2*steps;
    end
end
