function results = dualAxis(rho1,rho2,omega1,omega2,phi1_0,phi2_0,t,n)
% function results = dualAxis(rho1,rho2,omega1,omega2,phi1_0,phi2_0,t,n)
% Solves a dual-axis spiral problem.
% 
% Previous definitions:
%   - C: pole of the primary axis.
%   - S: pole of the secondary axis.
%   - P: sub-satellite point.
%   - E: pole of the Euler axis.
% Inputs:
%   - rho1: angular distance from C to S [deg].
%   - rho2: angular distance from S to P [deg].
%   - omega1: angular velocity of S about C [rad/s].
%   - omega2: angular velocity of P about C [rad/s].
%   - phi1_0: initial orientation of S relative to C [deg].
%   - phi2_0: azimuthal initial orientation of P relative to C [deg].
%   - t: time of interest [s], can be a vector of size 1xm.
%   - n (optional): number of parameters to be returned. The mininum is 1
%     and the maximum is 11. The default value is 11.
% Output:
%   - results: a matrix of size mxn containing the different computed
%     parameters for the specified times. For each row the order is:
%       1) Azimuth of S about C relative to the prime meridian [deg].
%       2) Azimuth of P about S relative to C [deg].
%       3) Elevation of P relative to C [deg].
%       4) Change in azimuth of P about C [deg].
%       5) Azimuth of P about C [deg].
%       6) Angle from C to E [deg].
%       7) Angle from P to E [deg].
%       8) Rate of rotation about E [rad/s].
%       9) Velocity of P [rad/s].
%      10) Change in direction of motion of P [deg].
%      11) Direction of motion of P [deg].
% 
% Author: Aleix Pinardell
% Version: 1.1.1
% Date: 29 February 2016
%

if nargin < 8
    n = 11;
end
results = zeros(length(t),n);
if n == 0; return; end

phi1 = phi1_0 + rad2deg(omega1*t);
results(:,1) = phi1;
if n == 1; return; end

phi2 = phi2_0 + rad2deg(omega2*t);
results(:,2) = phi2;
if n == 2; return; end

deltap = acosd(cosd(rho1).*cosd(rho2) + sind(rho1).*sind(rho2).*cosd(phi2));
delta = 90 - deltap;
results(:,3) = delta;
if n == 3; return; end

argument = (cosd(rho2)-cosd(rho1).*sind(delta))./(sind(rho1).*cosd(delta));
Dalpha = acos2d(validate(argument,-1,1), -Hd(phi2));
results(:,4) = Dalpha;
if n == 4; return; end

alpha = mod(phi1 + Dalpha, 360);
results(:,5) = alpha;
if n == 5; return; end

deltaEp = mod(atan2d(omega2.*sind(rho1),omega1+omega2.*cosd(rho1)), 180);
results(:,6) = deltaEp;
if n == 6; return; end

rhoE = acosd(sind(delta).*cosd(deltaEp) + cosd(delta).*sind(deltaEp).*cosd(Dalpha));
results(:,7) = rhoE;
if n == 7; return; end

omegaE = sqrt(omega1.^2 + omega2.^2 + 2*omega1.*omega2.*cosd(rho1));
results(:,8) = omegaE;
if n == 8; return; end

v = omegaE*sind(rhoE);
results(:,9) = v;
if n == 9; return; end

argument = (cosd(deltaEp)-cosd(rhoE).*sind(delta))./(sind(rhoE).*cosd(delta));
Dpsi = acos2d(validate(argument,-1,1), Hd(Dalpha));
results(:,10) = Dpsi;
if n == 10; return; end

psi = mod(Dpsi - 90, 360);
results(:,11) = psi;
