function [X,t] = nint(X0,f,tf,h,method,returnLast)
% function [X,t] = nint(X0,@dXdt,tf,h,method,returnLast)
% Integrates a state vector using the a numerical integration method.
% Inputs:
%   X0: initial state. Scalar or column vector.
%   f: name of the MATLAB function that evaluates the time derivative of
%      the state. This function must be of the form dXdt = f(t,X) and must
%      return a scalar/vector of the same size as X.
%   tf: the final time. Initial time is assumed to be zero.
%   h: the step-size integration time.
%   method: a string indicating the integration method to be used.
%           Currently supported: 'euler' and 'rk4'.
%   returnLast: (optional) 1 to force the function to return only the
%               evaluations at tf, 0 or unspecified to return all the evaluations.
% Outputs:
%   X: integrated state(s) (vector(s)).
%   t: integration time(s).
%
%  Author: Aleix Pinardell
%  Version: 1.2.1
%  Date: 6 December 2015

if nargin < 6 % if returnLast not specified, set it to 0
    returnLast = 0;
end
t = 0:h:tf;
if t(end) < tf % if last element of t is not tf, add it
    t = [t tf];
end
X = X0*ones(size(t)); % initialise the integrated state vector(s)
for i = 2:length(t)
    dt = (t(i) - t(i-1)); % = h in most cases, but maybe not for last step
    switch lower(method) % string in lowercase, if user writes 'Euler' also works
        case 'euler'
            Phi = f(t(i-1), X(:,i-1));
        case 'rk4'
            k1 = f(t(i-1), X(:,i-1));
            k2 = f(t(i-1) + dt/2, X(:,i-1) + dt*k1/2);
            k3 = f(t(i-1) + dt/2, X(:,i-1) + dt*k2/2);
            k4 = f(t(i-1) + dt, X(:,i-1) + dt*k3);
            Phi = 1/6*(k1 + 2*k2 + 2*k3 + k4);
        otherwise
            error('The integration method you specified is not supported.');
    end
    X(:,i) = X(:,i-1) + dt*Phi;
end
if returnLast
    t = t(end);
    X = X(:,end);
end
