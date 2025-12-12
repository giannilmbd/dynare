% Wrapper class for the PARDISO MEX files.
% In particular, ensures that there is no memory leak, thanks to the destructor.

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

classdef pardiso < handle

properties
    pt
    iparm
    dparm
end

methods
    function obj = pardiso(number_of_threads)
        if exist('pardiso_init') ~= 3 || exist('pardiso_free') ~= 3 || exist('pardiso_solve') ~= 3
            error('Cannot find the PARDISO MEX files. If you are compiling from source, pass -Dpardiso=enabled to meson and make sure that Panua PARDISO library is in the library search path.')
        end

        if ispc && matlab_ver_less_than('9.14')
            % The implementation of Intel OpenMP shipped with MATLAB < R2023a is too old to be
            % compatible with the one against which PARDISO library version 20240630 is linked
            warning('Your MATLAB version may be too old to be compatible with the PARDISO library provided by Panua. If you encounter crashes, you should upgrade to MATLAB R2023a at least.')
        end

        [obj.pt, obj.iparm, obj.dparm] = pardiso_init(number_of_threads);
    end

    function delete(obj)
        if ~isempty(obj.pt)
            pardiso_free(obj.pt, obj.iparm, obj.dparm)
            obj.pt = [];
        end
    end

    % Use standard solver
    function x = solve(obj, A, b)
        obj.iparm(4) = 0;
        [x, obj.pt, obj.iparm, obj.dparm] = pardiso_solve(A, b, obj.pt, obj.iparm, obj.dparm);
    end

    % Use CGS solver. Assumes that the standard solver has already been called on similar matrices.
    function x = solve_cgs(obj, A, b, iter_tol)
        obj.iparm(4) = -10*floor(log10(iter_tol)) + 1;
        [x, obj.pt, obj.iparm, obj.dparm] = pardiso_solve(A, b, obj.pt, obj.iparm, obj.dparm, 33);
    end
end

end
