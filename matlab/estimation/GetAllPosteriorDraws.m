function draws = GetAllPosteriorDraws(options_, dname, fname, column, FirstLine , FirstMhFile, blck)
% function draws = GetAllPosteriorDraws(options_, dname, fname, column, FirstLine , FirstMhFile, blck)
% Gets all posterior draws.
%
% INPUTS
% - options_               [struct]             Dynare's options.
% - dname                  [char]               name of directory with results.
% - fname                  [char]               name of mod file.
% - column                 [integer or string]  scalar, parameter index
%                                               'all' if all parameters are 
%                                               requested
% - FirstMhFile            [integer]            optional scalar, first MH file to be used
% - FirstLine              [integer]            optional scalar, first line in first MH file to be used
% - blck:                  [integer]  scalar, desired block to read (optional)
%
% OUTPUTS
% - draws:                 [double]   NumberOfDraws×1 vector, draws from posterior distribution.

% Copyright © 2005-2025 Dynare Team
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

if isnumeric(column) && column<0
    error('GetAllPosteriorDraws:: column cannot be negative');
end

if issmc(options_)
    if ishssmc(options_)
        % Load draws from the posterior distribution
        pfiles = dir(sprintf('%s/hssmc/particles-*.mat', dname));
        posterior = load(sprintf('%s/hssmc/particles-%u-%u.mat', dname, length(pfiles), length(pfiles)));
        if column>0 || strcmp(column,'all')
            if strcmp(column,'all')
                draws = transpose(posterior.particles);
            else
                draws = transpose(posterior.particles(column,:));
            end
        else
            draws = posterior.tlogpostkernel;
        end
    elseif isdsmh(options_)
        % Load draws from the posterior distribution
        posterior = load(sprintf('%s/dsmh/parameters_particles_final.mat', dname));
        if column>0 || strcmp(column,'all')
            if strcmp(column,'all')
                draws = transpose(posterior.particles);
            else
                draws = transpose(posterior.particles(column,:));
            end
        else
            draws = posterior.tlogpostkernel;
        end
    elseif isonline(options_)
        % Load draws from the posterior distribution
        posterior = load(sprintf('%s/online/parameters_particles_final.mat', dname));
        if column>0 || strcmp(column,'all')
            if strcmp(column,'all')
                draws = transpose(posterior.param);
            else
                draws = transpose(posterior.param(column,:));
            end
        else
            draws=NaN(size(posterior.param,2),1);
        end
    elseif isdime(options_)
        posterior = load(sprintf('%s%s%s%schains.mat', dname, filesep(), 'dime', filesep()));
        tune = posterior.tune;
        chains = posterior.chains(end-tune:end,:,:);
        if column>0 || strcmp(column,'all')
            chains = reshape(chains, [], size(chains, 3));
            if strcmp(column,'all')
                draws = chains;
            else
                draws = chains(:,column);
            end
        else
            draws = posterior.lprobs;
        end
    else
        error('GetAllPosteriorDraws:: case should not happen. Please contact the developers')
    end
else
    DirectoryName = CheckPath('metropolis',dname);
    record=load_last_mh_history_file(DirectoryName,fname);
    if nargin<5 || isempty(FirstLine)
        FirstLine = record.KeepedDraws.FirstLine;
    end
    if nargin<6 || isempty(FirstMhFile)
        FirstMhFile = record.KeepedDraws.FirstMhFile;
    end
    TotalNumberOfMhFile = sum(record.MhDraws(:,2));
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    [nblcks, npar] = size(record.LastParameters);
    iline = FirstLine;
    linee = 1;
    if nargin==7
        blocks_to_load=blck;
        nblcks=length(blck);
    else
        blocks_to_load=1:nblcks;
    end
    if strcmp(column,'all')
        draws = zeros(NumberOfDraws*nblcks,npar);
    else
        draws = zeros(NumberOfDraws*nblcks,1);
    end
    iline0=iline;
    if column>0
        for blck = blocks_to_load
            iline=iline0;
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
                NumberOfLines = size(x2(iline:end,:),1);
                if strcmp(column,'all')
                    draws(linee:linee+NumberOfLines-1,:) = x2(iline:end,:);
                else
                    draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
                end
                linee = linee+NumberOfLines;
                iline = 1;
            end
        end
    else
        for blck = blocks_to_load
            iline=iline0;
            for file = FirstMhFile:TotalNumberOfMhFile
                load([DirectoryName '/'  fname '_mh' int2str(file) '_blck' int2str(blck)],'logpo2')
                NumberOfLines = size(logpo2(iline:end),1);
                draws(linee:linee+NumberOfLines-1) = logpo2(iline:end);
                linee = linee+NumberOfLines;
                iline = 1;
            end
        end
    end
end
