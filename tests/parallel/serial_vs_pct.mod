/*
 * This file compares results from serial and parallel (PCT) Metropolis-Hastings estimations.
 * It checks if the x2 and logpo2 fields in the generated mat files are identical
 * between serial and parallel runs for different random seeds that were set
 * prior to the estimation command.
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

@#include "folder/endovars.mod" // model is irrelevant
@#include "_ls2003.inc"          // as we are running only MATLAB code below

@#define TESTS = ["start", "123", "default", "noseed"]

if ~isoctave && matlab.internal.parallel.isPCTInstalled

fprintf('\nTESTS BETWEEN SERIAL AND PCT ESTIMATIONS DONE IN OTHER MOD FILES:\n');

@#for t in TESTS
serial_1 = load(['serial_mhblocks_1_@{t}/metropolis/serial_mhblocks_1_mh1_blck1.mat']);
pct_1 = load(['pct_mhblocks_1_@{t}/metropolis/pct_mhblocks_1_mh1_blck1.mat']);

if isequaln(serial_1.x2, pct_1.x2) && isequaln(serial_1.logpo2, pct_1.logpo2)
    fprintf('✅ single blocks in serial and pct are equal for seed set to @{t}!\n');
else
    error('❌ single blocks in serial and pct are not equal for seed set to @{t}!')
end
@#endfor

@#for t in TESTS
serial_1 = load(['serial_mhblocks_2_@{t}/metropolis/serial_mhblocks_2_mh1_blck1.mat']);
serial_2 = load(['serial_mhblocks_2_@{t}/metropolis/serial_mhblocks_2_mh1_blck2.mat']);
pct_1 = load(['pct_mhblocks_2_@{t}/metropolis/pct_mhblocks_2_mh1_blck1.mat']);
pct_2 = load(['pct_mhblocks_2_@{t}/metropolis/pct_mhblocks_2_mh1_blck2.mat']);

if isequaln(serial_1.x2, pct_1.x2) && isequaln(serial_1.logpo2, pct_1.logpo2)
    fprintf('✅ first blocks in serial and pct are equal for seed set to @{t}!\n');
else
    error('❌ first blocks in serial and pct are not equal for seed set to @{t}!')
end
if isequaln(serial_2.x2, pct_2.x2) && isequaln(serial_2.logpo2, pct_2.logpo2)
    fprintf('✅ second blocks in serial and pct are equal for seed set to @{t}!\n');
else
    error('❌ second blocks in serial and pct are not equal for seed set to @{t}!')
end
@#endfor
fprintf('\n');

end