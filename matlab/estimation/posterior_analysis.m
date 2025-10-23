function oo_ = posterior_analysis(type,arg1,arg2,arg3,options_,M_,oo_,estim_params_)
% oo_ = posterior_analysis(type,arg1,arg2,arg3,options_,M_,oo_,estim_params_)
% Inputs
% - type            [string]        type of object to be computed
% - arg1                            first input argument of called function
%                                   (usually variable list)
% - arg2                            second input argument of called function
%                                   (usually variable or shock list)
% - arg3                            first input argument of called function
%                                   (nar or FEVD steps)
% - options_        [structure]     Dynare structure defining global options.
% - M_              [structure]     Dynare structure describing the model.
% - oo_             [structure]     Dynare structure where the results are saved.
% - estim_params_   [structure]     structure storing information about estimated
%                   parameters
% Outputs:
% - oo_             [structure]     Dynare structure where the results are saved.

% Copyright © 2008-2025 Dynare Team
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

info = check_posterior_analysis_data(type,M_,options_);
SampleSize = options_.sub_draws;
switch info
    case 0
        disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
        error('Check the options of the estimation command...')
    case {1,2} %need to get draws and store posterior
        MaxMegaBytes = options_.MaximumNumberOfMegaBytes;
        drsize = size_of_the_reduced_form_model(oo_.dr);
        if drsize*SampleSize>MaxMegaBytes
            drsize=0;
        end
        select_posterior_draws(M_,options_,oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state,estim_params_,SampleSize,drsize); %save draws to disk
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
    case {4,5} %process draws and save files
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
    case 6 %just read out
        [ivar,vartan] = get_variables_list(options_,M_);
        nvar = length(ivar);
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan);
    otherwise
        error('posterior_analysis:: Check_posterior_analysis_data gave a meaningless output!')
end



function oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan)
narg1 = 8;
narg2 = 10;
if ~(nargin==narg1 || nargin==narg2)
    error('posterior_analysis:: Call to function job is buggy!')
end
switch type
  case 'variance'
    if nargin==narg1
        [nvar,vartan] = ...
            dsge_simulated_theoretical_covariance(SampleSize,arg3,M_,options_,oo_,'posterior');
    end
    oo_ = covariance_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                 vartan,nvar,arg1,arg2,options_.mh_conf_sig,oo_,options_);
  case 'decomposition'
    if nargin==narg1
        [~,vartan] = ...
            dsge_simulated_theoretical_variance_decomposition(SampleSize,M_,options_,oo_,'posterior');
    end
    oo_ = variance_decomposition_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                             M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
    if ~all(diag(M_.H)==0)
        if strmatch(arg1,options_.varobs,'exact')
            observable_name_requested_vars=intersect(vartan,options_.varobs,'stable');
            oo_ = variance_decomposition_ME_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                                        [M_.exo_names;'ME'],arg2,observable_name_requested_vars,arg1,options_.mh_conf_sig,oo_,options_);
        end
    end
  case 'correlation'
    oo_ = correlation_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                  vartan,nvar,arg1,arg2,arg3,options_.mh_conf_sig,oo_,M_,options_);
  case 'conditional decomposition'
    if nargin==narg1
        [~,vartan] = ...
            dsge_simulated_theoretical_conditional_variance_decomposition(SampleSize,arg3,M_,options_,oo_,'posterior');
    end
    oo_ = conditional_variance_decomposition_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                                      arg3,M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
    if ~all(diag(M_.H)==0)
        if strmatch(arg1,options_.varobs,'exact')
            oo_ = conditional_variance_decomposition_ME_mc_analysis(SampleSize,'posterior',M_.dname,M_.fname,...
                                                              arg3,[M_.exo_names;'ME'],arg2,vartan,arg1,options_.mh_conf_sig,oo_,options_);
        end
    end
  otherwise
    disp('Not yet implemented')
end

end% job end



function SampleAddress = select_posterior_draws(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,estim_params_,SampleSize,drsize)
% SampleAddress = select_posterior_draws(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,estim_params_,SampleSize,drsize)
% --------------------------------------------------------------------------
% Selects a sample of draws from the posterior distribution and if nargin>1
% saves the draws in _pdraws mat files (metropolis folder). If drsize>0
% the dr structure, associated to the parameters, is also saved in _pdraws.
% This routine assures an _mh file cannot be opened twice.
%
% INPUTS
%   o M_                    [structure]     MATLAB's structure describing the model
%   o options_              [structure]     MATLAB's structure describing the current options
%   o dr                    [structure]     Reduced form model.
%   o endo_steady_state     [vector]        steady state value for endogenous variables
%   o exo_steady_state      [vector]        steady state value for exogenous variables
%   o exo_det_steady_state  [vector]        steady state value for exogenous deterministic variables                                    
%   o SampleSize            [integer]       Size of the sample to build.
%   o drsize                [double]        structure dr is drsize megaoctets.
%
% OUTPUTS
%   o SampleAddress  [integer]  A (SampleSize*4) array, each line specifies the
%                               location of a posterior draw:
%                                  Column 2 --> Chain number
%                                  Column 3 --> (mh) File number
%                                  Column 4 --> (mh) line number

% Number of parameters:
npar = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np;

% Select one task:
switch nargin
  case 8
    info = 0;
  case 9
    MAX_mega_bytes = 100;% Should be an option...
    if drsize>0
        info=2;
    else
        info=1;
    end
    drawsize = drsize+npar*8/1048576;
  otherwise
    error('select_posterior_draws:: Unexpected number of input arguments!')
end

if ~issmc(options_)
    MetropolisFolder = CheckPath('metropolis',M_.dname);
else
    if ishssmc(options_)
        [MetropolisFolder] = CheckPath('hssmc',M_.dname);
    elseif isdime(options_)
        [MetropolisFolder] = CheckPath('dime',M_.dname);
    else
        error('check_posterior_analysis_data:: case should not happen. Please contact the developers')
    end
end

ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

parameter_draws=get_posterior_subsample(M_,options_,SampleSize);

% Selected draws in the posterior distribution, and if drsize>0
% reduced form solutions, are saved on disk.
if info
    %delete old stale files before creating new ones
    delete_stale_file([BaseName '_posterior_draws*.mat'])
    NumberOfDrawsPerFile = fix(MAX_mega_bytes/drawsize);
    NumberOfFiles = ceil(SampleSize*drawsize/MAX_mega_bytes);
    NumberOfLinesLastFile = SampleSize - (NumberOfFiles-1)*NumberOfDrawsPerFile;
    linee = 0;
    fnum  = 1;
    if fnum < NumberOfFiles
        pdraws = cell(NumberOfDrawsPerFile,info);
    else
        pdraws = cell(NumberOfLinesLastFile,info);
    end
    for i=1:SampleSize
        linee = linee+1;
        pdraws(i,1) = {parameter_draws(i,:)};
        if info==2
            M_ = set_parameters_locally(M_,pdraws{linee,1});
            [dr,~,M_.params] = compute_decision_rules(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            pdraws(linee,2) = { dr };
        end
        if (fnum < NumberOfFiles && linee == NumberOfDrawsPerFile) || (fnum <= NumberOfFiles && i==SampleSize)
            linee = 0;
            save([BaseName '_posterior_draws' num2str(fnum) '.mat'],'pdraws','estim_params_')
            fnum = fnum+1;
            if fnum < NumberOfFiles
                pdraws = cell(NumberOfDrawsPerFile,info);
            else
                pdraws = cell(NumberOfLinesLastFile,info);
            end
        end
    end
end

end % select_posterior_draws end

end % posterior_analysis end
