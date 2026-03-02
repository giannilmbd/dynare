/*
 * Runs a serial estimation with 2 Metropolis-Hastings blocks.
 */

/*
 * Copyright © 2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

@#include "folder/endovars.mod"
@#include "_ls2003.inc"

@#define MH_NBLOCKS = 2
@#define DIRNAME    = "serial_mhblocks_" + ((string) MH_NBLOCKS)

if ~isoctave && matlab.internal.parallel.isPCTInstalled
options_.parallel_info.use_pct.estimation.sampler = false;
delete(gcp('nocreate')); % make sure no parallel pool is open
parpool('local', 2); % open pool to check whether desire of user to not use PCT is honored and pool is restored correctly on cleanup
pool_before = gcp('nocreate');
pool_before_NumWorkers = pool_before.NumWorkers;
pool_before_Profile = pool_before.Cluster.Profile;
@#include "_estimation_commands_with_different_seeds.inc"

% Test that parallel pool was correctly restored after estimation
pool_after = gcp('nocreate');
if isempty(pool_after)
    error('❌ Parallel pool was not restored after serial estimation!');
end
if pool_after.NumWorkers ~= pool_before_NumWorkers
    error('❌ Parallel pool was restored with wrong number of workers: expected %d, got %d', pool_before_NumWorkers, pool_after.NumWorkers);
end
if ~strcmp(pool_after.Cluster.Profile, pool_before_Profile)
    error('❌ Parallel pool was restored with wrong profile: expected ''%s'', got ''%s''', pool_before_Profile, pool_after.Cluster.Profile);
end
fprintf('✅ Parallel pool correctly restored after serial estimation (%d workers, %s profile)\n', pool_after.NumWorkers, pool_after.Cluster.Profile);
delete(pool_after); % clean up
end
