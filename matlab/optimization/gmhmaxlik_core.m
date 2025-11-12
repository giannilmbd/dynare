function [PostMod,PostVar,Scale,PostMean,fcount] = gmhmaxlik_core(ObjFun,xparam1,mh_bounds,options,iScale,info,MeanPar,VarCov,varargin)
% [PostMod,PostVar,Scale,PostMean,fcount] = gmhmaxlik_core(ObjFun,xparam1,mh_bounds,options,iScale,info,MeanPar,VarCov,varargin)
% (Dirty) Global minimization routine of (minus) a likelihood (or posterior density) function.
%
% INPUTS
%   o ObjFun     [char]     string specifying the name of the objective function.
%   o xparam1    [double]   (p*1) vector of parameters to be estimated.
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters.
%   o options    [structure] options for the optimization algorithm (options_.gmhmaxlik).
%   o iScale     [double]   scalar specifying the initial of the jumping distribution's scale parameter.
%   o info       [char]     string, empty or equal to 'LastCall'.
%   o MeanPar    [double]   (p*1) vector specifying the initial posterior mean.
%   o VarCov     [double]   (p*p) matrix specifying the initial posterior covariance matrix.
%   o gend       [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data       [double]   (T*n) matrix of data ==> varargin{2}.
%
% OUTPUTS
%   o PostMod    [double]   (p*1) vector, evaluation of the posterior mode.
%   o PostVar    [double]   (p*p) matrix, evaluation of the posterior covariance matrix.
%   o Scale      [double]   scalar specifying the scale parameter that should be used in
%                           an eventual Metropolis-Hastings algorithm.
%   o PostMean   [double]   (p*1) vector, evaluation of the posterior mean.
%   o fcount     [integer]  scalar, number of function evaluations.
%
% ALGORITHM
%   Metropolis-Hastings with an constantly updated covariance matrix for
%   the jump distribution. The posterior mean, variance and mode are
%   updated (in step 2) with the following rules:
%
%   \[
%       \mu_t = \mu_{t-1} + \frac{1}{t}\left(\theta_t-\mu_{t-1}\right)
%   \]
%
%   \[
%       \Sigma_t = \Sigma_{t-1} + \mu_{t-1}\mu_{t-1}'-\mu_{t}\mu_{t}' +
%                  \frac{1}{t}\left(\theta_t\theta_t'-\Sigma_{t-1}-\mu_{t-1}\mu_{t-1}'\right)
%   \]
%
%   and
%
%   \[
%       \mathrm{mode}_t = \left\{
%                       \begin{array}{ll}
%                         \theta_t, & \hbox{if } p(\theta_t|\mathcal Y) > p(\mathrm{mode}_{t-1}|\mathcal Y) \\
%                         \mathrm{mode}_{t-1}, & \hbox{otherwise.}
%                       \end{array}
%                     \right.
%   \]
%
%   where $t$ is the iteration, $\mu_t$ the estimate of the posterior mean
%   after $t$ iterations, $\Sigma_t$ the estimate of the posterior
%   covariance matrix after $t$ iterations, $\mathrm{mode}_t$ is the
%   evaluation of the posterior mode after $t$ iterations and
%   $p(\theta_t|\mathcal Y)$ is the posterior density of parameters
%   (specified by the user supplied function "fun").
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

fcount = 0;
npar = length(xparam1);

NumberOfIterations = options.number;
MaxNumberOfTuningSimulations   = options.nscale;
MaxNumberOfClimbingSimulations = options.nclimb;
AcceptanceTarget               = options.target;
console_mode                   = options.console_mode;
CovJump = VarCov;
ModePar = xparam1;

%% [1] I tune the scale parameter.
if ~options.silent
    [hh_fig, length_of_old_string] = wait_bar.run(0,[],'Tuning of the scale parameter...', console_mode, 0, 'Tuning of the scale parameter.');
end
j = 1; jj  = 1;
isux = 0; jsux = 0; test = 0;
ix2 = ModePar;% initial condition!
ilogpo2 = - feval(ObjFun,ix2,varargin{:});% initial posterior density
fcount = fcount + 1;
mlogpo2 = ilogpo2;
try
    dd = transpose(chol(CovJump));
catch
    dd = eye(length(CovJump));
end
while j<=MaxNumberOfTuningSimulations
    proposal = iScale*dd*randn(npar,1) + ix2;
    if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
        logpo2 = - feval(ObjFun,proposal,varargin{:});
        fcount = fcount + 1;
    else
        logpo2 = -inf;
    end
    % I move if the proposal is enough likely...
    if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
        ix2 = proposal;
        if logpo2 > mlogpo2
            ModePar = proposal;
            mlogpo2 = logpo2;
        end
        ilogpo2 = logpo2;
        isux = isux + 1;
        jsux = jsux + 1;
    end % ... otherwise I don't move.
    prtfrc = j/MaxNumberOfTuningSimulations;
    if ~options.silent && mod(j, 10)==0
        [~,length_of_old_string]=wait_bar.run(prtfrc,hh_fig,sprintf('Acceptance ratio [during last 500]: %f [%f]',isux/j,jsux/jj),console_mode,length_of_old_string);
    end
    if  j/500 == round(j/500)
        test1 = jsux/jj;
        cfactor = test1/AcceptanceTarget;
        if cfactor>0
            iScale = iScale*cfactor;
        else
            iScale = iScale/10;
        end
        jsux = 0; jj = 0;
        if cfactor>0.90 && cfactor<1.10
            test = test+1;
        end
        if test>4
            break
        end
    end
    j = j+1;
    jj = jj + 1;
end

if ~options.silent
    wait_bar.close(hh_fig,console_mode);
end
%% [2] One block metropolis, I update the covariance matrix of the jumping distribution
if ~options.silent
    [hh_fig, length_of_old_string] = wait_bar.run(0,[],'Metropolis-Hastings...', console_mode, 0, 'Estimation of the posterior covariance...');
end
j = 1;
isux = 0;
ilogpo2 = - feval(ObjFun,ix2,varargin{:});
fcount = fcount + 1;
while j<= NumberOfIterations
    j = j+1;
    proposal = iScale*dd*randn(npar,1) + ix2;
    if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
        logpo2 = - feval(ObjFun,proposal,varargin{:});
        fcount = fcount + 1;
    else
        logpo2 = -inf;
    end
    % I move if the proposal is enough likely...
    if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
        ix2 = proposal;
        if logpo2 > mlogpo2
            ModePar = proposal;
            mlogpo2 = logpo2;
        end
        ilogpo2 = logpo2;
        isux = isux + 1;
        jsux = jsux + 1;
    end % ... otherwise I don't move.
    prtfrc = j/NumberOfIterations;
    if ~options.silent && mod(j, 10)==0
        [~,length_of_old_string]=wait_bar.run(prtfrc,hh_fig,sprintf('Acceptance ratio: %f',isux/j),console_mode,length_of_old_string);
    end
    % I update the covariance matrix and the mean:
    oldMeanPar = MeanPar;
    MeanPar = oldMeanPar + (1/j)*(ix2-oldMeanPar);
    CovJump = CovJump + oldMeanPar*oldMeanPar' - MeanPar*MeanPar' + ...
              (1/j)*(ix2*ix2' - CovJump - oldMeanPar*oldMeanPar');
end
if ~options.silent
    wait_bar.close(hh_fig,console_mode);
end
PostVar = CovJump;
PostMean = MeanPar;
%% [3 & 4] I tune the scale parameter (with the new covariance matrix) if
%% this is the last call to the routine, and I climb the hill (without
%% updating the covariance matrix)...
if strcmpi(info,'LastCall')
    if ~options.silent
        [hh_fig, length_of_old_string] = wait_bar.run(0,[],'Tuning of the scale parameter...', console_mode, 0, 'Tuning of the scale parameter.');
    end
    j = 1; jj  = 1;
    isux = 0; jsux = 0;
    test = 0;
    ilogpo2 = - feval(ObjFun,ix2,varargin{:});% initial posterior density
    fcount = fcount + 1;
    dd = transpose(chol(CovJump));
    while j<=MaxNumberOfTuningSimulations
        proposal = iScale*dd*randn(npar,1) + ix2;
        if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
            logpo2 = - feval(ObjFun,proposal,varargin{:});
            fcount = fcount + 1;
        else
            logpo2 = -inf;
        end
        % I move if the proposal is enough likely...
        if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
            ix2 = proposal;
            if logpo2 > mlogpo2
                ModePar = proposal;
                mlogpo2 = logpo2;
            end
            ilogpo2 = logpo2;
            isux = isux + 1;
            jsux = jsux + 1;
        end % ... otherwise I don't move.
        prtfrc = j/MaxNumberOfTuningSimulations;
        if ~options.silent && mod(j, 10)==0
            [~,length_of_old_string]=wait_bar.run(prtfrc,hh_fig,sprintf('Acceptance ratio [during last 1000]: %f [%f]',isux/j,jsux/jj),console_mode,length_of_old_string);
        end
        if j/1000 == round(j/1000)
            test1 = jsux/jj;
            cfactor = test1/AcceptanceTarget;
            iScale = iScale*cfactor;
            jsux = 0; jj = 0;
            if cfactor>0.90 && cfactor<1.10
                test = test+1;
            end
            if test>4
                break
            end
        end
        j = j+1;
        jj = jj + 1;
    end
    if ~options.silent
        wait_bar.close(hh_fig,console_mode);
    end
    Scale = iScale;
    %%
    %% Now I climb the hill
    %%
    if options.nclimb
        if ~options.silent
            [hh_fig, length_of_old_string] = wait_bar.run(0,[],' ', console_mode, 0, 'Now I am climbing the hill...');
        end
        j = 1; jj  = 1;
        jsux = 0;
        test = 0;
        while j<=MaxNumberOfClimbingSimulations
            proposal = iScale*dd*randn(npar,1) + ModePar;
            if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
                logpo2 = - feval(ObjFun,proposal,varargin{:});
                fcount = fcount + 1;
            else
                logpo2 = -inf;
            end
            if logpo2 > mlogpo2% I move if the proposal is higher...
                ModePar = proposal;
                mlogpo2 = logpo2;
                jsux = jsux + 1;
            end % otherwise I don't move...
            prtfrc = j/MaxNumberOfClimbingSimulations;
            if ~options.silent && mod(j, 10)==0
                [~,length_of_old_string]=wait_bar.run(prtfrc,hh_fig,sprintf('%f Jumps / MaxStepSize %f',jsux,sqrt(max(diag(iScale*CovJump)))),console_mode,length_of_old_string);
            end
            if  j/200 == round(j/200)
                if jsux<=1
                    test = test+1;
                else
                    test = 0;
                end
                jsux = 0;
                jj = 0;
                if test>4% If I do not progress enough I reduce the scale parameter
                         % of the jumping distribution (cooling down the system).
                    iScale = iScale/1.10;
                end
                if sqrt(max(diag(iScale*CovJump)))<10^(-4)
                    break% Steps are too small!
                end
            end
            j = j+1;
            jj = jj + 1;
        end
        if ~options.silent
            wait_bar.close(hh_fig,console_mode);
        end
    end %climb
else
    Scale = iScale;
end
PostMod = ModePar;