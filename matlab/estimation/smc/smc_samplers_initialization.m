function [particles, tlogpostkernel, loglikelihood] = smc_samplers_initialization(objective_function, sampler, n, Prior, SimulationFolder, nsteps, base_seed, parallel_info)
% function [particles, tlogpostkernel, loglikelihood] = smc_samplers_initialization(objective_function, sampler, n, Prior, SimulationFolder, nsteps, base_seed, parallel_info)
% Initialize SMC samplers by drawing initial particles in the prior distribution.
%
% INPUTS
% - objective_function   [char]     string specifying the name of the objective function (posterior kernel).
% - sampler              [char]     name of the sampler.
% - n                    [integer]  scalar, number of particles.
% - Prior                [class]    prior information
% - SimulationFolder     [char]     name of output folder
% - nsteps               [integer]  number of steps
% - base_seed            [integer]  base seed for the random number generator (options_.DynareRandomStreams.seed)
% - parallel_info        [struct]   parallel info structure (options_.parallel_info)
%
% OUTPUTS
% - particles             [double]    p×n matrix of particles
% - tlogpostkernel        [double]    n×1 vector of posterior kernel values for the particles
% - loglikelihood         [double]    n×1 vector of likelihood values for the particles
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2022-2026 Dynare Team
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

% SMC does not support Dragonfly (parallel command line option)
if parallel_info.isHybridMatlabOctave
    error('smc_samplers_initialization:dragonfly_not_supported', ...
        'SMC samplers do not support the parallel command line option with mixed MATLAB/Octave workers.')
end

dprintf('Estimation:%s: Initialization...', sampler)

% Delete old mat files storing chains and particles if any...
chains_file = sprintf('%s%schains.mat', SimulationFolder, filesep());
particle_files = sprintf('%s%sparticles*.mat', SimulationFolder, filesep());
if isfile(chains_file) || ~isempty(dir(particle_files))
    if isfile(chains_file); delete(chains_file); end
    if ~isempty(dir(particle_files)); delete(particle_files); end
    dprintf('Estimation:%s: Old %s-files successfully erased.', sampler, sampler)
end

dprintf('Estimation:%s: Searching for initial values...', sampler);
particles = zeros(Prior.length(), n);
tlogpostkernel = zeros(n, 1);
loglikelihood = zeros(n, 1);

% Pre-compute per-particle RNG seeds so that serial and parallel produce identical results (same as in posterior_sampler_initialization)
InitialSeeds = struct();
for j=1:n
    set_dynare_seed_local_options([], false, base_seed+j);
    [InitialSeeds(j).Unifor, InitialSeeds(j).Normal, InitialSeeds(j).global_stream] = get_dynare_random_generator_state();
end

% Setup parallel execution using the shared utility function setup_parallel_execution
% which handles PCT availability check, user preferences, pool management, and cleanup.
[run_with_pct, restore_pool] = setup_parallel_execution(parallel_info.use_pct.estimation.smc_initialization, 'smc_samplers_initialization'); %#ok<ASGLU>
% restore_pool holds an onCleanup object that restores the pool state when this function ends or crashes; we need to keep it in the scope of this function and so ignore the warning about unused variable

% Simulate a pool of particles characterizing the prior distribution (with the additional constraint that the likelihood is finite)
t0 = tic;
if run_with_pct
    parfor j=1:n
        set_dynare_random_generator_state(InitialSeeds(j).Unifor, InitialSeeds(j).Normal, InitialSeeds(j).global_stream);
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
        set_dynare_random_generator_state(InitialSeeds(j).Unifor, InitialSeeds(j).Normal, InitialSeeds(j).global_stream);
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
% Restore RNG state on the main process so that callers see a consistent
% state regardless of whether parfor or for was used.
set_dynare_random_generator_state(InitialSeeds(n).Unifor, InitialSeeds(n).Normal, InitialSeeds(n).global_stream);
tt = toc(t0);
save(sprintf('%s%sparticles-1-%u.mat', SimulationFolder, filesep(), nsteps), 'particles', 'tlogpostkernel', 'loglikelihood')
dprintf('Estimation:%s: Initial values found (%.2fs)', sampler, tt)
skipline()
