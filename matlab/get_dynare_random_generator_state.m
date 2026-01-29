function [state_u,state_n,global_stream] = get_dynare_random_generator_state()
% [state_u,state_n,global_stream] = get_dynare_random_generator_state()
% Get state of MATLAB/Octave random generator 
% Outputs:
%   - state_u           [integer]       state of uniform RNG
%   - state_n           [integer]       state of normal RNG
%   - global_stream     [randstream]    current default RNG stream
%
% Notes: 
%  - Octave, like older versions of MATLAB keeps one generator for uniformly distributed numbers and
%    one for normally distributed numbers. 
%  - For compatibility, we return two vectors, which are identical for MATLAB. 
%  - For MATLAB, we also return the current stream as it may change when
%       using parfor

% Copyright © 2010-2026 Dynare Team
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

if ~isoctave
    global_stream = RandStream.getGlobalStream();
    state_u = global_stream.State;
    state_n = state_u;
else
    global_stream = [];
    state_u = rand('state');
    state_n = randn('state');
end