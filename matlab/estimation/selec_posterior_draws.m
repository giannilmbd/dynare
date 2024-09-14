function SampleAddress = selec_posterior_draws(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,estim_params_,SampleSize,drsize)
% Selects a sample of draws from the posterior distribution and if nargin>1
% saves the draws in _pdraws mat files (metropolis folder). If drsize>0
% the dr structure, associated to the parameters, is also saved in _pdraws.
% This routine assures an _mh file cannot be opened twice.
%
% INPUTS
%   o M_                    [structure]     Matlab's structure describing the model
%   o options_              [structure]     Matlab's structure describing the current options
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
%
% SPECIAL REQUIREMENTS
%   None.
%

% Copyright © 2006-2024 Dynare Team
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

% Number of parameters:
npar = estim_params_.nvx;
npar = npar + estim_params_.nvn;
npar = npar + estim_params_.ncx;
npar = npar + estim_params_.ncn;
npar = npar + estim_params_.np;

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
    error('selec_posterior_draws:: Unexpected number of input arguments!')
end

if ~issmc(options_)
    MetropolisFolder = CheckPath('metropolis',M_.dname);
else
    if ishssmc(options_)
        [MetropolisFolder] = CheckPath('metropolis',M_.dname);
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