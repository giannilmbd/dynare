function set_dynare_random_generator_state(state_u,state_n,current_stream)
%Set state of MATLAB/Octave random generator 
% Outputs:
%   - None
% Inputs:
%   - state_u           [integer]       state of uniform RNG
%   - state_n           [integer]       state of normal RNG
%   - global_stream     [randstream]    RNG stream to be made default stream
%
% Notes: 
%  - Octave, like older versions of MATLAB keeps one generator for uniformly distributed numbers and
%    one for normally distributed numbers. 
%  - For compatibility, we input two vectors, which need to be identical for MATLAB. 
%  - For MATLAB, we also meed to define the current stream as it may have changed when
%       parfor was used

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
    if nargin<3
        current_stream = RandStream.getGlobalStream();
    end
    if ~isequal(state_u,state_n)
        error(['You are using the new MATLAB RandStream mechanism ' ...
            'with a single random generator, but the values ' ...
            'of the state of the uniformly ' ...
            'distributed numbers and of the state of the ' ...
            'normally distributed numbers are different. Something must be ' ...
            'wrong, such as reloading old Metropolis runs, ' ...
            'computed on a different version of MATLAB. If you ' ...
            'don''t understand the origin of the problem, ' ...
            'please, contact Dynare''s development team.'])
    end
    current_stream.State = state_u;
    RandStream.setGlobalStream(current_stream);
else
    rand('state',state_u);
    randn('state',state_n);
end
