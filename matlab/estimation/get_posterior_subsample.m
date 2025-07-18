function [pdraws, log_posterior]=get_posterior_subsample(M_,options_,SampleSize)
% function [pdraws, log_posterior]=get_posterior_subsample(M_,options_,SampleSize)
% Checks the status of posterior analysis and in particular if files need to be
% created or updated; called by posterior_analysis.m
%
% Inputs:
%   M_              [structure]     Dynare model structure
%   options_        [structure]     Dynare options structure
%   SampleSize      [integer]       number of desired subdraws
%
% Outputs:
%   pdraws          [double]        SampleSize by npar matrix of parameter draws
%   log_posterior   [double]        SampleSize by 1 vector of log posterior

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

if ~issmc(options_)
    % Get information about the MCMC:
    record=load_last_mh_history_file([M_.dname filesep 'metropolis'], M_.fname);
    npar=size(record.InitialParameters,2);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine;
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

    % Randomly select draws in the posterior distribution:
    SampleAddress = zeros(SampleSize,4);
    for iter = 1:SampleSize
        ChainNumber = ceil(rand*options_.mh_nblck);
        DrawNumber  = ceil(rand*NumberOfDraws);
        SampleAddress(iter,1) = DrawNumber;
        SampleAddress(iter,2) = ChainNumber;
        if DrawNumber <= MAX_nruns-FirstLine+1 % in first file where FirstLine may be bigger than 1
            MhFileNumber = FirstMhFile;
            MhLineNumber = FirstLine+DrawNumber-1;
        else
            DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
            MhFileNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns);
            MhLineNumber = DrawNumber-(MhFileNumber-FirstMhFile-1)*MAX_nruns;
        end
        SampleAddress(iter,3) = MhFileNumber;
        SampleAddress(iter,4) = MhLineNumber;
    end
    SampleAddress = sortrows(SampleAddress,[3 2]);

    % Selected draws in the posterior distribution, and if drsize>0
    % reduced form solutions, are saved on disk.
    pdraws = NaN(SampleSize,npar);
    if nargout>1
        log_posterior = NaN(SampleSize,1);
    end
    old_mhfile = 0;
    old_mhblck = 0;
    for iter = 1:SampleSize
        mhfile = SampleAddress(iter,3);
        mhblck = SampleAddress(iter,2);
        if (mhfile ~= old_mhfile) || (mhblck ~= old_mhblck)
            temp=load([M_.dname filesep 'metropolis' filesep M_.fname '_mh' num2str(mhfile) '_blck' num2str(mhblck) '.mat'],'x2','logpo2');
        end
        pdraws(iter,:) = temp.x2(SampleAddress(iter,4),:);
        if nargout>1
            log_posterior(iter,1) = temp.logpo2(SampleAddress(iter,4));
        end
        old_mhfile = mhfile;
        old_mhblck = mhblck;
    end
else
    if ishssmc(options_)
        % Load draws from the posterior distribution
        pfiles = dir(sprintf('%s/hssmc/particles-*.mat', M_.dname));
        posterior = load(sprintf('%s/hssmc/particles-%u-%u.mat', M_.dname, length(pfiles), length(pfiles)));
        NumberOfDraws=size(posterior.particles,2);
        DrawNumbers  = ceil(rand(SampleSize,1)*NumberOfDraws);
        pdraws = transpose(posterior.particles(:,DrawNumbers));
        if nargout>1
            log_posterior = posterior.tlogpostkernel(DrawNumbers);
        end
    elseif isdime(options_)
        posterior = load(sprintf('%s%s%s%schains.mat', M_.dname, filesep(), 'dime', filesep()));
        tune = posterior.tune;
        chains = reshape(posterior.chains(end-tune:end,:,:), [], size(posterior.chains, 3));
        NumberOfDraws=size(chains,1);
        DrawNumbers  = ceil(rand(SampleSize,1)*NumberOfDraws);
        pdraws = chains(DrawNumbers,:);
        if nargout>1
            log_posterior = posterior.lprobs(DrawNumbers);
        end
    end
end