function [PostMode, HessianMatrix, Scale, ModeValue, fcount] = gmhmaxlik(fun, xinit, Hinit, iscale, bounds, priorstd, gmhmaxlikOptions, OptimizationOptions, varargin)
% [PostMode, HessianMatrix, Scale, ModeValue, fcount] = gmhmaxlik(fun, xinit, Hinit, iscale, bounds, priorstd, gmhmaxlikOptions, OptimizationOptions, varargin)
% Prepare (Dirty) Global minimization routine of (minus) a likelihood (or posterior density) function
%
% INPUTS
%   o fun                 [char]      string specifying the name of the objective function.
%   o xinit               [double]    (p*1) vector of parameters to be estimated (initial values).
%   o Hinit               [double]    (p*p) matrix specifying the initial Hessian matrix (if empty, the prior covariance matrix is used).
%   o iscale              [double]    scalar specifying the initial of the jumping distribution's scale parameter.
%   o bounds              [double]    (p*2) matrix defining lower and upper bounds for the parameters.
%   o priorstd            [double]    (p*1) vector specifying the standard deviations of the prior (if empty, the prior covariance matrix is used).
%   o gmhmaxlikOptions    [structure]  options for the optimization algorithm (options_.gmhmaxlik).
%   o OptimizationOptions [cell array] user-defined options for the optimization algorithm.
%   o varargin            [cell array] additional arguments for the objective function.
%
% OUTPUTS
%   o PostMode      [double]   (p*1) vector, evaluation of the posterior mode.
%   o HessianMatrix [double]   (p*p) matrix, evaluation of the Hessian matrix at the posterior mode.
%   o Scale         [double]   scalar specifying the scale parameter that should be used in an eventual Metropolis-Hastings algorithm.
%   o ModeValue     [double]   scalar, function value at the posterior mode.
%   o fcount        [integer]  scalar, number of function evaluations.
%
% ALGORITHM
%  See gmhmaxlik_core.m for details.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2006-2025 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% Set default options

fcount = 0;
if ~isempty(Hinit)
    gmhmaxlikOptions.varinit = 'previous';
else
    gmhmaxlikOptions.varinit = 'prior';
end

if ~isempty(OptimizationOptions)
    DynareOptionslist = read_key_value_string(OptimizationOptions);
    for i=1:rows(DynareOptionslist)
        switch DynareOptionslist{i,1}
          case 'NumberOfMh'
            gmhmaxlikOptions.iterations = DynareOptionslist{i,2};
          case 'ncov-mh'
            gmhmaxlikOptions.number = DynareOptionslist{i,2};
          case 'nscale-mh'
            gmhmaxlikOptions.nscale = DynareOptionslist{i,2};
          case 'nclimb-mh'
            gmhmaxlikOptions.nclimb = DynareOptionslist{i,2};
          case 'InitialCovarianceMatrix'
            switch DynareOptionslist{i,2}
              case 'previous'
                if isempty(Hinit)
                    error('gmhmaxlik: No previous estimate of the Hessian matrix available! You cannot use the InitialCovarianceMatrix option!')
                else
                    gmhmaxlikOptions.varinit = 'previous';
                end
              case {'prior', 'identity'}
                gmhmaxlikOptions.varinit = DynareOptionslist{i,2};
              otherwise
                error('gmhmaxlik: Unknown value for option ''InitialCovarianceMatrix''!')
            end
          case 'AcceptanceRateTarget'
            gmhmaxlikOptions.target = DynareOptionslist{i,2};
            if gmhmaxlikOptions.target>1 || gmhmaxlikOptions.target<eps
                error('gmhmaxlik: The value of option AcceptanceRateTarget should be a double between 0 and 1!')
            end
          otherwise
            warning(['gmhmaxlik: Unknown option (' DynareOptionslist{i,1}  ')!'])
        end
    end
end

% Evaluate the objective function.
OldModeValue = feval(fun,xinit,varargin{:});
fcount = fcount + 1;

if ~exist('MeanPar','var')
    MeanPar = xinit;
end

switch gmhmaxlikOptions.varinit
  case 'previous'
    CovJump = inv(Hinit);
  case 'prior'
    % The covariance matrix is initialized with the prior
    % covariance (a diagonal matrix) %%Except for infinite variances ;-)
    stdev = priorstd;
    indx = find(isinf(stdev));
    stdev(indx) = ones(length(indx),1)*sqrt(10);
    vars = stdev.^2;
    CovJump = diag(vars);
  case 'identity'
    vars = ones(length(priorstd),1)*0.1;
    CovJump = diag(vars);
  otherwise
    error('gmhmaxlik: This is a bug! Please contact the developers.')
end

OldPostVariance = CovJump;
OldPostMean = xinit;
OldPostMode = xinit;
Scale = iscale;

for i=1:gmhmaxlikOptions.iterations
    if i<gmhmaxlikOptions.iterations
        flag = '';
    else
        flag = 'LastCall';
    end
    [PostMode, PostVariance, Scale, PostMean, fc] = gmhmaxlik_core(fun, OldPostMode, bounds, gmhmaxlikOptions, Scale, flag, MeanPar, OldPostVariance, varargin{:});
    fcount = fcount + fc;
    ModeValue = feval(fun, PostMode, varargin{:});
    fcount = fcount + 1;
    dVariance = max(max(abs(PostVariance-OldPostVariance)));
    dMean = max(abs(PostMean-OldPostMean));
    if ~gmhmaxlikOptions.silent
        skipline()
        printline(58,'=')
        disp(['   Change in the posterior covariance matrix = ' num2str(dVariance) '.'])
        disp(['   Change in the posterior mean = ' num2str(dMean) '.'])
        disp(['   Current mode = ' num2str(ModeValue)])
        disp(['   Mode improvement = ' num2str(abs(OldModeValue-ModeValue))])
        disp(['   New value of jscale = ' num2str(Scale)])
        printline(58,'=')
    end
    OldModeValue = ModeValue;
    OldPostMean = PostMean;
    OldPostVariance = PostVariance;
end

HessianMatrix = inv(PostVariance);

if ~gmhmaxlikOptions.silent
    skipline()
    disp(['Optimal value of the scale parameter = ' num2str(Scale)])
    skipline()
end