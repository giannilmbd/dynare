function C = pagemtimes(A, B)
% PAGEMTIMES  Page-wise matrix multiplication (Octave compatibility shim)
%
% C = pagemtimes(A, B) performs page-wise matrix multiplication of 3-D arrays.
% For pages i, C(:,:,i) = A(:,:,i) * B(:,:,i)
%
% This is a compatibility function for Octave. MATLAB R2020b+ has this built-in.

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

    % Get dimensions
    [m, n, p] = size(A);
    [n2, k, p2] = size(B);

    if n ~= n2
        error('pagemtimes:DimensionMismatch', ...
              'Inner dimensions must agree for matrix multiplication');
    end

    if p ~= p2 && p ~= 1 && p2 ~= 1
        error('pagemtimes:PageMismatch', ...
              'Number of pages must match or one array must have a single page');
    end

    % Handle singleton expansion
    if p == 1 && p2 > 1
        A = repmat(A, [1, 1, p2]);
        p = p2;
    elseif p2 == 1 && p > 1
        B = repmat(B, [1, 1, p]);
    end

    % Preallocate output
    C = zeros(m, k, p);

    % Perform page-wise multiplication
    for i = 1:p
        C(:, :, i) = A(:, :, i) * B(:, :, i);
    end
end
