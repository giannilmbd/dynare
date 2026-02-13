/*
 * Runs a serial estimation with 1 Metropolis-Hastings block.
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

@#define MH_NBLOCKS = 1
@#define DIRNAME    = "serial_mhblocks_" + ((string) MH_NBLOCKS)

if ~isoctave && matlab.internal.parallel.isPCTInstalled
options_.parallel_info.use_pct.estimation = false;
delete(gcp('nocreate')); % make sure no parallel pool is open
@#include "_estimation_commands_with_different_seeds.inc"
end
