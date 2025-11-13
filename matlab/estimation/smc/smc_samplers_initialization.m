function [particles, tlogpostkernel, loglikelihood] = smc_samplers_initialization(objective_function, sampler, n, Prior, SimulationFolder, nsteps)
% function [particles, tlogpostkernel, loglikelihood] = smc_samplers_initialization(objective_function, sampler, n, Prior, SimulationFolder, nsteps)
% Initialize SMC samplers by drawing initial particles in the prior distribution.
%
% INPUTS
% - objective_function [char]     string specifying the name of the objective function (posterior kernel).
% - sampler          [char]     name of the sampler.
% - n                [integer]  scalar, number of particles.
% - Prior            [class]    prior information
% - SimulationFolder [char]     name of output folder
% - nsteps           [integer]  number of steps
%
% OUTPUTS
% - particles             [double]    p×n matrix of particles
% - tlogpostkernel        [double]    n×1 vector of posterior kernel values for the particles
% - loglikelihood         [double]    n×1 vector of likelihood values for the particles
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2022-2025 Dynare Team
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

dprintf('Estimation:%s: Initialization...', sampler)

% Delete old mat files storing particles if any...
matfiles = sprintf('%s%sparticles*.mat', SimulationFolder, filesep());
files = dir(matfiles);
if ~isempty(files)
    delete(matfiles);
    dprintf('Estimation:%s: Old %s-files successfully erased.', sampler, sampler)
end

dprintf('Estimation:%s: Searching for initial values...', sampler);
particles = zeros(Prior.length(), n);
tlogpostkernel = zeros(n, 1);
loglikelihood = zeros(n, 1);

% Simulate a pool of particles characterizing the prior distribution (with the additional constraint that the likelihood is finite)
set_dynare_seed('default');

t0 = tic;
if ~isoctave && matlab.internal.parallel.isPCTInstalled
    sc = parallel.pool.Constant(RandStream('Threefry'));
    parfor j=1:n
        stream = sc.Value;
        stream.Substream = j;
        RandStream.setGlobalStream(stream); % set the seed in each iteration of parfor loops
        notvalid = true;
        while notvalid
            candidate = Prior.draw();
            if Prior.admissible(candidate)
                particles(:,j) = candidate;
                [tlogpostkernel(j), loglikelihood(j)] = tempered_likelihood(objective_function, candidate, 0.0, Prior);
                if isfinite(loglikelihood(j)) % if returned log-density is Inf or Nan (penalized value)
                    notvalid = false;
                end
            end
        end
    end
else
    for j=1:n
        notvalid = true;
        while notvalid
            candidate = Prior.draw();
            if Prior.admissible(candidate)
                particles(:,j) = candidate;
                [tlogpostkernel(j), loglikelihood(j)] = tempered_likelihood(objective_function, candidate, 0.0, Prior);
                if isfinite(loglikelihood(j)) % if returned log-density is Inf or Nan (penalized value)
                    notvalid = false;
                end
            end
        end
    end
end
tt = toc(t0);
save(sprintf('%s%sparticles-1-%u.mat', SimulationFolder, filesep(), nsteps), 'particles', 'tlogpostkernel', 'loglikelihood')
dprintf('Estimation:%s: Initial values found (%.2fs)', sampler, tt)
skipline()
