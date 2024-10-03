function [sub_draws, error_flag]=set_number_of_subdraws(M_,options_)
% function [sub_draws, error_flag]=set_number_of_subdraws(M_,options_)
% Set option field sub_draws based on size of available sample
% Inputs:
%   M_          [structure]     Dynare model structure
%   options_    [structure]     Dynare options structure
%
% Outputs:
%   sub_draws   [scalar]        number of sub-draws
%   error_flag  [boolean]       error indicate of insufficient draws are
%                               available

% Copyright © 2024 Dynare Team
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


error_flag=false;
if ~issmc(options_)
    record=load_last_mh_history_file([M_.dname filesep 'metropolis'], M_.fname);
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDrawsPerChain = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    NumberOfDraws=NumberOfDrawsPerChain*record.Nblck;
    if isempty(options_.sub_draws)
        sub_draws = min(options_.posterior_max_subsample_draws, ceil(NumberOfDraws));
    else
        if options_.sub_draws>NumberOfDraws*record.Nblck
            skipline()
            disp(['The value of option sub_draws (' num2str(options_.sub_draws) ') is greater than the number of available draws in the MCMC (' num2str(NumberOfDraws*options_.mh_nblck) ')!'])
            disp('You can either change the value of sub_draws, reduce the value of mh_drop, or run another mcmc (with the load_mh_file option).')
            skipline()
            error_flag = true;
            sub_draws=options_.sub_draws; %pass through
        else
            sub_draws=options_.sub_draws; %pass through
        end
    end
else
    if ishssmc(options_)
        % Load draws from the posterior distribution
        pfiles = dir(sprintf('%s/hssmc/particles-*.mat', M_.dname));
        posterior = load(sprintf('%s/hssmc/particles-%u-%u.mat', M_.dname, length(pfiles), length(pfiles)));
        NumberOfDraws = size(posterior.particles,2);
    elseif isdime(options_)
        posterior = load(sprintf('%s%s%s%schains.mat', M_.dname, filesep(), 'dime', filesep()));
        tune = posterior.tune;
        chains = reshape(posterior.chains(end-tune:end,:,:), [], size(posterior.chains, 3));
        NumberOfDraws=size(chains,1);
    else
        error('set_number_of_subdraws:: case should not happen. Please contact the developers')
    end
    if isempty(options_.sub_draws)
        sub_draws = min(options_.posterior_max_subsample_draws, NumberOfDraws);
    else
        if options_.sub_draws>NumberOfDraws
            skipline()
            disp(['The value of option sub_draws (' num2str(options_.sub_draws) ') is greater than the number of available draws in the MCMC (' num2str(NumberOfDraws) ')!'])
            disp('You can either change the value of sub_draws or run another MCMC.')
            skipline()
            error_flag = true;
            sub_draws=options_.sub_draws; %pass through
        else
            sub_draws=options_.sub_draws; %pass through
        end
    end
end