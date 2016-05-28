function [x,f,eval] = newrap(x0,F,dF,tol)
% function [x,f,eval] = newrap(x0,F,dF,tol)
% Finds the minimum of a function F(x) using the Newton-Raphson method.
% Inputs:
%   x0: initial value.
%   F: the function to be evalutaed. Can be a MATLAB function handle of the
%      type F = @(x) ... or the name of a MATLAB function with the
%      following syntax: f = F(x).
%   dF: the derivative of F. Can also be a function handle or a function.
%   tol: the iteration process stops when f differs less than tol from zero.
% Outputs:
%   x: the root of F.
%   f: the value of F at x (close to zero).
%   eval: the total number of F and dF evaluations.
%
%  Author: Aleix Pinardell
%  Version: 1.1
%  Date: 12 December 2015

x = x0;
eval = 0;
f = 2*tol;
while abs(f) >= tol
    f = F(x);
    df = dF(x);
    x = x - f/df;
    eval = eval + 2;
end
