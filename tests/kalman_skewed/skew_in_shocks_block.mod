% Tests for preprocessor capabilities of skew in (stochastic) shocks block
% Note: manual tests for error messages by uncommenting specific lines
% -------------------------------------------------------------------------

% Copyright © 2025 Dynare Team
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

@#include "ireland2004_model.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 0: check different ways to specify skewness and co-skewness %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
% Single shock notation
skew eta_a, eta_a, eta_a = 1;
% skew eta_a = 1; % same as above
skew eta_e = 2;
skew eta_r = 3;
skew eta_z = 4;
%skew eta_u = 5; % triggers error because of undeclared variable (eta_u)
%skew yhat = 6; % triggers error because skewness only allowed for exogenous
%skew ALPHA_PI = 7; % triggers error because skewness only allowed for exogenous

% Two same shocks, one different
skew eta_a, eta_e, eta_a = 10;
skew eta_e, eta_r, eta_e = 11;
skew eta_z, eta_a, eta_z = 12;
skew eta_r, eta_r, eta_z = 13;
%skew eta_r, eta_r, eta_u = 14; % triggers error because of undeclared variable (eta_u) [nostrict works]
%skew eta_r, eta_r, yhat = 15; % triggers error because skewness only allowed for exogenous
%skew eta_r, ALPHA_PI, eta_e = 16; % triggers error because skewness only allowed for exogenous

% All three shocks different
skew eta_a, eta_e, eta_r = 20;
%skew eta_e, eta_a, eta_r = 20; % Same as above, just different order
%skew eta_r, eta_a, eta_e = 20; % Same as above, just different order
skew eta_a, eta_e, eta_z = 21;
skew eta_a, eta_r, eta_z = 22;
skew eta_e, eta_r, eta_z = 23;

% Fill remaining two-same-one-different families
skew eta_a, eta_a, eta_z = 31;
skew eta_a, eta_a, eta_r = 32;
skew eta_e, eta_e, eta_a = 33;
skew eta_e, eta_e, eta_z = 34;
skew eta_z, eta_z, eta_e = 35;
skew eta_z, eta_z, eta_r = 36;
skew eta_r, eta_r, eta_a = 37;
% Intentionally leave one family unset so its permutations remain zero
% skew eta_r, eta_r, eta_e = 38
end;

fprintf('\n======= M_.Skew_e =======\n');
disp(M_.Skew_e);
fprintf('\n=========================\n');

fprintf('\n=== Testing Co-Skewness Matrix Permutations ===\n\n');
fprintf('Matrix dimensions: %d x %d x %d\n\n', size(M_.Skew_e));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: Check symmetry across all permutations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 1: Checking permutation symmetry...\n');
errors = 0;
total_checks = 0;
    
for i = 1:size(M_.Skew_e, 1)
    for j = 1:size(M_.Skew_e, 1)
        for k = 1:size(M_.Skew_e, 1)
            % Get all 6 permutations
            val_ijk = M_.Skew_e(i,j,k);
            val_ikj = M_.Skew_e(i,k,j);
            val_jik = M_.Skew_e(j,i,k);
            val_jki = M_.Skew_e(j,k,i);
            val_kij = M_.Skew_e(k,i,j);
            val_kji = M_.Skew_e(k,j,i);
            
            % Collect all permutations
            perms = [val_ijk, val_ikj, val_jik, val_jki, val_kij, val_kji];
            
            % Check if all are equal (within numerical tolerance)
            if ~iszero(perms - perms(1))
                errors = errors + 1;
                fprintf('  ERROR: Permutations not equal for (%s, %s, %s)\n', ...
                    M_.exo_names{i}, M_.exo_names{j}, M_.exo_names{k});
                fprintf('    (%d,%d,%d) = %.6f\n', i, j, k, val_ijk);
                fprintf('    (%d,%d,%d) = %.6f\n', i, k, j, val_ikj);
                fprintf('    (%d,%d,%d) = %.6f\n', j, i, k, val_jik);
                fprintf('    (%d,%d,%d) = %.6f\n', j, k, i, val_jki);
                fprintf('    (%d,%d,%d) = %.6f\n', k, i, j, val_kij);
                fprintf('    (%d,%d,%d) = %.6f\n', k, j, i, val_kji);
            end
            total_checks = total_checks + 1;
        end
    end
end
    
if errors == 0
    fprintf('  ✓ PASSED: All %d permutations are consistent!\n\n', total_checks);
else
    error('  ✗ FAILED: Found %d inconsistencies out of %d checks\n\n', errors, total_checks);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: Display non-zero unique elements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 2: Summary of non-zero skewness values...\n');
reported = zeros(size(M_.Skew_e, 1), size(M_.Skew_e, 1), size(M_.Skew_e, 1));
count = 0;

for i = 1:size(M_.Skew_e, 1)
    for j = 1:size(M_.Skew_e, 1)
        for k = 1:size(M_.Skew_e, 1)
            if abs(M_.Skew_e(i,j,k)) > 0 && reported(i,j,k) == 0
                count = count + 1;
                % Mark all permutations as reported
                reported(i,j,k) = 1;
                reported(i,k,j) = 1;
                reported(j,i,k) = 1;
                reported(j,k,i) = 1;
                reported(k,i,j) = 1;
                reported(k,j,i) = 1;
                
                fprintf('  [%d] Skew(%s, %s, %s) = %.6f\n', ...
                    count, M_.exo_names{i}, M_.exo_names{j}, M_.exo_names{k}, M_.Skew_e(i,j,k));
            end
        end
    end
end

if count == 0
    fprintf('  No non-zero elements found.\n');
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: Check diagonal elements (own-shock skewness) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 3: Own-shock skewness (diagonal elements)...\n');
for i = 1:size(M_.Skew_e, 1)
    val = M_.Skew_e(i,i,i);
    if abs(val) > 0
        fprintf('  Skew(%s, %s, %s) = %.6f\n', ...
            M_.exo_names{i}, M_.exo_names{i}, M_.exo_names{i}, val);
    end
end
fprintf('\n');

fprintf('=== Tests Completed ===\n\n');