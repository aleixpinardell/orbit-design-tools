function [x_best,f_best] = geneticAlgorithm(f,range,varargin)
% function [x_best,f_best] = geneticAlgorithm(f,range,varargin)
% Finds the optimum of a function f(x) using a genetic algorithm.
% Inputs:
%   f: objective function. Can be either a MATLAB function handle of the
%      type f =@(x) ... or the name of a MATLAB function with the syntax
%      value = f(x). The input parameter x and the output of the function
%      f are row-vectors of the same size, 1xm.
%   range: a matrix of size 2xm containing the lower and upper limits of
%          the different parameters, which defines the search area.
%   varargin: a cell array containing pairs of keys and values, used to
%             define the behaviour of the genetic algorithm. The allowed
%             keys and their corresponding allowed values are:
%             - 'Objective': 'min' or 'max'. Default is 'min'.
%             - 'PopulationSize': default is 100.
%             - 'BitResolution': row-vector of integers with each column
%                corresponding to the bit resolution of a parameter, or a
%                integer, in which case the same bit resolution is used for
%                all the parameters. Default is 10.
%             - 'Elitism': percentage of parents that are copied directly
%                to the next generation. Default is 0.1.
%             - 'ParentSelection': 'random' or 'fittest'. Default is
%                'random'.
%             - 'MutationProbability': default is 1e-4.
%             - 'ConvergenceThreshold': maximum change in value of the
%                best-fitted solution of the current generation with
%                respect to that of the previous one for which two
%                consecutive generations are said to yield similar results.
%                Default is 1e-3.
%             - 'PlateauCount': minimum number of consecutive generations
%                with similar results for which the process can be
%                considered to have reached plateau. Default is 10.
%             - 'GenerationLimit': maximum number of generations. Default
%                is 100.
%             - 'History': 1 (true) to return values of all generations or
%                0 (false) to return only those of last generation.
% Outputs:
%   x_best: (a history of) the values of the parameters that render f optimum.
%   f_best: (a history of) the optimum value of f.
% 
% Author: Aleix Pinardell
% Version: 1.1
% Date: 7 February 2015


% Check valid size of range
[nvar,two] = size(range);
if two ~= 2
    error('The size for range parameter is not valid.');
end

% Define valid arguments
okargs = {'Objective', 'PopulationSize', 'BitResolution', 'Elitism', ...
    'ParentSelection', 'MutationProbability', 'ConvergenceThreshold', ...
    'PlateauCount', 'GenerationLimit', 'History'};

% Set default values
maximize = 0;
popsize = 100;
bitres = 10*ones(nvar,1);
fittest = 0.1;
selectmode = 'random';
mutprob = 1e-4;
tol = 1e-3;
genconv = 10;
maxit = 100;
history = 0;

% Read input arguments
for j=1:2:(nargin-2)
    pname = varargin{j};
    pval = varargin{j+1};
    pname = validatestring(pname,okargs,'geneticAlgorithm');
    switch(pname)
        case 'Objective'
            pval = validatestring(pval,{'max','min'});
            if strcmp(pval,'max')
                maximize = 1;
            end
        case 'PopulationSize'
            if isnumeric(pval) && isequal(round(pval),pval) && pval >= 1
                popsize = pval;
            else
                error('Value of PopulationSize must be an integer larger or equal than 1.');
            end
        case 'Elitism'
            if isnumeric(pval) && pval >= 0 && pval < 1
                fittest = pval;
            else
                error('Value of Elitism must be a positive real number lower than 1.');
            end
        case 'BitResolution'
            if isnumeric(pval) && isequal(round(pval),pval)
                [nbit,mbit] = size(pval);
                if ( nbit == nvar && mbit == 1 ) || ( mbit == nvar && nbit == 1 )
                    if nbit == 1
                        pval = pval';
                    end
                    bitres = pval;
                elseif size(pval) == [1 1]
                    bitres = pval*ones(nvar,1);
                else
                    error('Unconsistent size for BitResolution parameter.');
                end
            else
                error('Value of BitResolution must be an integer (array).');
            end
        case 'ParentSelection'
            selectmode = validatestring(pval,{'random','fittest'});
        case 'MutationProbability'
            if isnumeric(pval) && pval >= 0 pval < 1
                mutprob = pval;
            else
                error('Value of MutationProbability must be a positive real number lower than 1.');
            end
        case 'ConvergenceThreshold'
            if isnumeric(pval)
                tol = pval;
            else
                error('Value of ConvergenceThreshold must be real positive number.');
            end
        case 'PlateauCount'
            if isnumeric(pval) && isequal(round(pval),pval) && pval >= 2
                genconv = pval;
            else
                error('Value of PlateauCount not valid (must be an integer larger or equal than 2).');
            end
        case 'GenerationLimit'
            if isnumeric(pval) && isequal(round(pval),pval)
                maxit = pval;
            else
                error('Value of GenerationLimit must be a positive integer.');
            end
        case 'History'
            if pval == 0 || pval == 1
                history = pval;
            else
                error('Value of History not valid (must be either 0 for false or 1 for ture).');
            end
    end
end

if maximize
    f = @(x) -f(x);
end

% Create initial random population
population = randi([0 1],[popsize sum(bitres)]);
vars = zeros(popsize,nvar);

nit = 0;
err = 1;
f_best = [];
x_best = [];
while nit < maxit
    % Determine values of variables for current population
    for i = 1:popsize
        k = 1;
        for j = 1:nvar
            b = bitres(j);
            kend = k + b - 1;
            bits = population(i,k:kend);
            vars(i,j) = range(j,1) + (range(j,2)-range(j,1))/(2^b-1)*2.^(0:b-1)*bits';
            k = kend + 1;
        end
    end
    
    % Determine values of function for current population
    fvals = f(vars);
    [f_best_new,index] = min(fvals);
    
    % Find best-fitted
    [~, indexes] = sort(fvals); % best-fitted first
    
    % Off-spring (crossover)
    if strcmp(selectmode,'random')
        ps = indexes(floor(fittest*popsize) + randperm(ceil((1-fittest)*popsize)));
    else
        ps = indexes(floor(fittest*popsize) + (1:ceil((1-fittest)*popsize)));
    end
    if length(ps)/2 ~= round(length(ps)/2) % check if number of parents is odd
        ps(1) = []; % fittest goes directly to next generation
    end
    for i = 1:2:length(ps)
        temp = population(ps(i),:);
        cindex = randi([2 length(temp)-1]);
        population(ps(i),cindex:end) = population(ps(i+1),cindex:end);
        population(ps(i+1),cindex:end) = temp(cindex:end);
    end
    
    % Prepare for next iteration
    nit = nit + 1;
    imprv = abs((f_best - f_best_new)/f_best_new);
    convergence = 0;
    for i = 1:length(imprv)
        if imprv(end+1-i) < 0.1
            convergence = convergence + 1;
        else
            break
        end
    end
    f_best = [f_best f_best_new];
    x_best = [x_best; vars(index,:)];
    if convergence >= genconv
        break
    end
    
    % Mutations
    for i = max(ceil(fittest*popsize),1):popsize
        parent = population(indexes(i),:);
        for j = 1:length(parent)
            if rand <= mutprob
                bit = population(indexes(i),j);
                population(indexes(i),j) = mod(bit+1,2); % 0->1 or 1->0
            end
        end
    end
end
if maximize
    f_best = -f_best;
end
if ~history
    f_best = f_best(end);
    x_best = x_best(end,:);
end
