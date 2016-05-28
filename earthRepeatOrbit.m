function parameter = earthRepeatOrbit(j,k,a,e,i,basic)
% function parameter = earthRepeatOrbit(j,k,a,e,i,basic)
% Considers a satellite in an Earth repeat orbit completing j orbits every
% k days. The numerical value of two of the following orbital parameters
% must be specified:
%   a (semi-major axis)
%   e (eccentricity)
%   i (inclination)
% The non-specified parameter must be introduced as an empty vector [].
% The function solves for the value of the non-specified parameter that
% satisfies the requirements on the other two parameters, j and k.
% 
% Inputs:
%   j: orbital revolutions required for the groundtrack to repeat [orbital revolutions]
%   k: sideral days required for the groundtrack to repeat [sideral days]
%   a: semi-major axis [km]
%   e: eccentricity [-]
%   i: inclination [deg]
%   basic: (optional) 1 to ignore the effects of J_2 on the argument of
%          pericenter and the mean anomaly, do not specify to consider them
% Output:
%   parameter: either a, e or i in the corresponding units
% 
% Author: Aleix Pinardell
% Version: 1.2
% Date: 27 November 2015
%

% Load pre-defined constants
constants;

if nargin < 6
    basic = 0;
end

% Determine which parameter so solve for (a: 0, e: 1, i: 2)
solveFor = 0;
if isempty(e)
    solveFor = 1;
elseif isempty(i)
    solveFor = 2;
end

% Find the number of different orbits to be considered and make all the
% input orbital parameters of consistent size.
n_orbits = max([length(j) length(k) length(a) length(e) length(i)]);
j = rowVectorOfSize(j,n_orbits);
k = rowVectorOfSize(k,n_orbits);
a = rowVectorOfSize(a,n_orbits);
e = rowVectorOfSize(e,n_orbits);
i = rowVectorOfSize(i,n_orbits);
i = deg2rad(i);

L_dot = 360;
L_2_ind_basic = -3*pi*J_2*R_E^2; % L_2 part independent of a, e and i
L_2_ind = -3*J_2*sqrt(mu_E)*R_E^2*180/pi*D_sideral; % same but for basic=0
if solveFor == 0 % Find semi-major axis
    % Estimate initial value of semi-major axis (J_2 ignored)
    a = ((D_sideral*k./(2*pi*j)).^2*mu_E).^(1/3);
    for o = 1:n_orbits
        a_corr = 1;
        a_tol = 1e-10;
        max_it = 100;
        it = 0;
        while a_corr > a_tol
            if basic % Consider effects of J_2 only on raan
                L_2 = L_2_ind_basic*cos(i(o))/(a(o)^2*(1-e(o)^2)^2);
                L_1 = -k(o)*2*pi/j(o) - L_2;
                a_new = ((-L_1*D_sideral/(2*pi)^2)^2*mu_E)^(1/3);
            else % Consider effects of J_2 on raan, omega and M
                common = L_2_ind*a(o)^-3.5;
                raan_dot = 0.5*common*cos(i(o))*(1-e(o)^2)^-2;
                omega_dot = -0.25*common*(5*cos(i(o))^2-1)*(1-e(o)^2)^-2;
                M_dot = -0.25*common*(3*cos(i(o))^2-1)*(1-e(o)^2)^-1.5;
                n = j(o)/k(o)*(L_dot - raan_dot) - (omega_dot + M_dot);
                a_new = (mu_E/(n/D_sideral*pi/180)^2)^(1/3);
            end
            a_corr = abs(a_new - a(o));
            a(o) = a_new;
            it = it + 1;
            if it > max_it
                fprintf('Solution did not converge (j=%g, k=%g, e=%g, i=%g).\n',...
                    j(o),k(o),e(o),rad2deg(i(o)));
                break
            end
        end
        % Uncomment these lines for debug purposes
        % raan_dot
        % omega_dot
        % M_dot
    end
    parameter = a;
else
    L_1 = -(2*pi)^2*sqrt(a.^3/mu_E)/D_sideral;
    L_2 = -k*2*pi./j - L_1;
    if solveFor == 1 % Find eccentricity
        if basic
            e = sqrt(1 - sqrt(L_2_ind_basic.*cos(i)./(a.^2.*L_2)));
        else
            % First e estimation: circular orbit
            e = rowVectorOfSize(0,n_orbits);
            for o = 1:n_orbits
                e_corr = 1;
                e_tol = 1e-10;
                max_it = 100;
                it = 0;
                while e_corr > e_tol
                    n = sqrt(mu_E/a(o)^3)/(pi/180/D_sideral);
                    common = L_2_ind*a(o)^-3.5;
                    omega_dot = -0.25*common*(5*cos(i(o))^2-1)*(1-e(o)^2)^-2;
                    M_dot = -0.25*common*(3*cos(i(o))^2-1)*(1-e(o)^2)^-1.5;
                    raan_dot = L_dot - (n + omega_dot + M_dot)*k(o)/j(o);
                    e_new = sqrt(1 - sqrt(0.5*common*cos(i(o))/raan_dot));
                    e_corr = abs(e_new - e(o));
                    e(o) = e_new;
                    it = it + 1;
                    if it > max_it
                        fprintf('Solution did not converge (j=%g, k=%g, a=%g, i=%g).\n',...
                            j(o),k(o),a(o),rad2deg(i(o)));
                        break
                    end
                end
            end
        end
        parameter = e;
    else % Find inclination
        if basic
            i = acos(L_2./L_2_ind_basic.*(1-e.^2).^2.*a.^2);
        else
            n = sqrt(mu_E./a.^3)/(pi/180/D_sideral);
            common = L_2_ind*a.^-3.5;
            c_omega = -0.25*common.*(1-e.^2).^-2;
            c_M = -0.25*common.*(1-e.^2).^-1.5;
            c_c = j./k*L_dot - n + c_omega + c_M; % term independent of cos(i)
            c_b = -j./k*0.5.*common.*(1-e.^2).^-2; % term in cos(i)
            c_a = -(5*c_omega + 3*c_M); % term in cos(i)^2
            cosi = (-c_b + sqrt(c_b.^2 - 4*c_a.*c_c))./(2*c_a);
            i = acos(cosi);
        end
        parameter = rad2deg(i);
    end
end
