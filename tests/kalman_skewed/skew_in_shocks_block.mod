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

% Helper to look up skewness value from sparse Skew_e (sorts indices, returns 0 if not found)
skew_lookup = @(S,i,j,k) sum(S(S(:,1)==min([i j k]) & S(:,2)==median([i j k]) & S(:,3)==max([i j k]), 4));

fprintf('\n======= M_.Skew_e (sparse 4-column format, sorted indices) =======\n');
disp(M_.Skew_e);
fprintf('\n=========================\n');

fprintf('\n=== Testing Co-Skewness Sparse Matrix ===\n\n');
fprintf('Number of stored entries: %d\n\n', size(M_.Skew_e, 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: Check that indices are stored in increasing order %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 1: Checking sorted index order...\n');
errors = 0;

for row = 1:size(M_.Skew_e, 1)
    i = M_.Skew_e(row, 1);
    j = M_.Skew_e(row, 2);
    k = M_.Skew_e(row, 3);

    if ~(i <= j && j <= k)
        errors = errors + 1;
        fprintf('  ERROR: Row %d has unsorted indices (%d, %d, %d)\n', row, i, j, k);
    end
end

if errors == 0
    fprintf('  PASSED: All %d entries have sorted indices (i <= j <= k)!\n\n', size(M_.Skew_e, 1));
else
    error('  FAILED: Found %d rows with unsorted indices\n\n', errors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: Check that any permutation lookup returns correct %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 2: Checking that any-permutation lookup works...\n');
errors = 0;
total_checks = 0;

for row = 1:size(M_.Skew_e, 1)
    i = M_.Skew_e(row, 1);
    j = M_.Skew_e(row, 2);
    k = M_.Skew_e(row, 3);
    v = M_.Skew_e(row, 4);

    % Look up all 6 permutations via skew_lookup (which sorts internally)
    perm_vals = [skew_lookup(M_.Skew_e,i,j,k), skew_lookup(M_.Skew_e,i,k,j), ...
                 skew_lookup(M_.Skew_e,j,i,k), skew_lookup(M_.Skew_e,j,k,i), ...
                 skew_lookup(M_.Skew_e,k,i,j), skew_lookup(M_.Skew_e,k,j,i)];

    if ~iszero(perm_vals - v)
        errors = errors + 1;
        fprintf('  ERROR: Lookup mismatch for (%s, %s, %s)\n', ...
            M_.exo_names{i}, M_.exo_names{j}, M_.exo_names{k});
    end
    total_checks = total_checks + 1;
end

if errors == 0
    fprintf('  PASSED: All %d entries return consistent values for any permutation!\n\n', total_checks);
else
    error('  FAILED: Found %d inconsistencies out of %d checks\n\n', errors, total_checks);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: Display non-zero unique elements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 3: Summary of non-zero skewness values...\n');

for count = 1:size(M_.Skew_e, 1)
    i = M_.Skew_e(count, 1);
    j = M_.Skew_e(count, 2);
    k = M_.Skew_e(count, 3);
    fprintf('  [%d] Skew(%s, %s, %s) = %.6f\n', ...
        count, M_.exo_names{i}, M_.exo_names{j}, M_.exo_names{k}, M_.Skew_e(count, 4));
end

if size(M_.Skew_e, 1) == 0
    fprintf('  No non-zero elements found.\n');
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 4: Check diagonal elements (own-shock skewness) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test 4: Own-shock skewness (diagonal elements)...\n');
for i = 1:M_.exo_nbr
    val = skew_lookup(M_.Skew_e, i, i, i);
    if abs(val) > 0
        fprintf('  Skew(%s, %s, %s) = %.6f\n', ...
            M_.exo_names{i}, M_.exo_names{i}, M_.exo_names{i}, val);
    end
end
fprintf('\n');

fprintf('=== Tests Completed ===\n\n');