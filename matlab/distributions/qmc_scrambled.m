function x = qmc_scrambled(nshocks, Nsam, scramble)
% x = qmc_scrambled(nshocks, Nsam, scramble)
%
% Generates scrambled Sobol quasi-Monte Carlo (QMC) sequences.
%
% INPUTS
% - nshocks                 [integer]   Number of dimensions (shocks).
% - Nsam                    [integer]   Number of samples to generate.
% - scramble                [integer]   Scrambling method:
%                                       0 = column-wise permutation (preserves QMC properties)
%                                       1 = column and row-wise permutation (preserves QMC properties)
%                                       2 = full scrambling (breaks QMC properties, equivalent to LHS)
%
% OUTPUTS
% - x                       [double]    Matrix of QMC draws (Nsam x nshocks).
%
% REMARKS
% Column-wise and row-wise scrambling preserve the quasi-Monte Carlo properties of the sequence,
% while element-wise scrambling breaks QMC properties but provides a Latin Hypercube Sampling equivalent.

% Copyright © 2025-2026 Dynare Team
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

lpmat = qmc_sequence(nshocks, int64(1), 0, Nsam)';
x=lpmat;
switch scramble
    case 0
        % column wise, keeps QMC properties
        icol = randperm(nshocks);
        x=lpmat(:,icol);
    case 1
        % column and row wise, keeps QMC properties
        iraw= randperm(Nsam);
        icol = randperm(nshocks);
        x=lpmat(:,icol);
        x=x(iraw,:);
    case 2
        % breaks QMC properties, equivalent to to lhs
        for k=1:nshocks
            iscramble = randperm(Nsam);
            x(:,k)=lpmat(iscramble,k);
        end
end