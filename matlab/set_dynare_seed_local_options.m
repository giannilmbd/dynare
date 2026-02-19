function DynareRandomStreams=set_dynare_seed_local_options(DynareRandomStreams,isHybridMatlabOctave,a,b)
% DynareRandomStreams=set_dynare_seed_local_options(DynareRandomStreams,isHybridMatlabOctave,a,b)
% Set seeds depending on MATLAB (octave) version
% Inputs:
%   o DynareRandomStreams   options structure relating to options_.DynareRandomStreams,
%                           potentially empty for initialization (a=='default')
%   o isHybridMatlabOctave  bool indicating whether parallel mixed pool is
%                           used
%   o a                     first input argument,
%                           for single argument input, either
%                               [number]       seed
%                               [string]        'default'
%                                               'reset'
%       	                                      'clock'
%                           for two inputs:
%                               [string]        indicating feasible RNG algorithm (default: 'mt19937ar')
%   o b                     second input argument
%                               [number]        seed for algorithm
%                                               specified with first input
% Outputs:
%   o DynareRandomStreams   options structure relating to options_.DynareRandomStreams


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

if nargin<3
    error('set_dynare_seed:: I need at least three input arguments!')
end

matlab_random_streams = ~(isoctave || isHybridMatlabOctave);

if matlab_random_streams% Use new MATLAB interface.
    if nargin==3
        if ischar(a) && strcmpi(a,'default')
            DynareRandomStreams.algo = 'mt19937ar';
            DynareRandomStreams.seed = 0;
            DynareRandomStreams.global_stream = RandStream(DynareRandomStreams.algo,'Seed',DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(DynareRandomStreams.global_stream));
            return
        end
        if ischar(a) && strcmpi(a,'reset')
            DynareRandomStreams.global_stream = RandStream(DynareRandomStreams.algo,'Seed',DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(DynareRandomStreams.global_stream));
            return
        end
        if ~ischar(a) || (ischar(a) && strcmpi(a, 'clock'))
            DynareRandomStreams.algo = 'mt19937ar';
            if ischar(a)
                DynareRandomStreams.seed = rem(floor(now*24*60*60), 2^32);
            else
                DynareRandomStreams.seed = a;
            end
            DynareRandomStreams.global_stream = RandStream(DynareRandomStreams.algo,'Seed',DynareRandomStreams.seed);
            reset(RandStream.setGlobalStream(DynareRandomStreams.global_stream));
            return
        end
        error('set_dynare_seed:: something is wrong in the calling sequence!')
    elseif nargin==4
        if ~ischar(a) || ~( strcmpi(a,'mcg16807') || ...
                            strcmpi(a,'mlfg6331_64') || ...
                            strcmpi(a,'mrg32k3a') || ...
                            strcmpi(a,'mt19937ar') || ...
                            strcmpi(a,'shr3cong') || ...
                            strcmpi(a,'swb2712') )
            disp('set_dynare_seed:: third argument must be string designing the uniform random number algorithm!')
            RandStream.list
            skipline()
            disp('set_dynare_seed:: Change the third input accordingly...')
            skipline()
            error(' ')
        end
        if ~isint(b)
            error('set_dynare_seed:: The fourth input argument must be an integer!')
        end
        DynareRandomStreams.algo = a;
        DynareRandomStreams.seed = b;
        DynareRandomStreams.global_stream = RandStream(DynareRandomStreams.algo,'Seed',DynareRandomStreams.seed);
        reset(RandStream.setGlobalStream(DynareRandomStreams.global_stream));
    end
else% Use old MATLAB interface.
    DynareRandomStreams.global_stream=[];
    if nargin==3
        if ischar(a) && strcmpi(a,'default')
            if isoctave
                DynareRandomStreams.algo = 'state';
            else
                DynareRandomStreams.algo = 'twister';
            end
            DynareRandomStreams.seed = 0;
            rand(DynareRandomStreams.algo,DynareRandomStreams.seed);
            randn('state',DynareRandomStreams.seed);
            return
        end
        if isempty(DynareRandomStreams) || ~isfield(DynareRandomStreams,'algo')
        %make sure algorithm is set
            if isoctave
                DynareRandomStreams.algo = 'state';
            else
                DynareRandomStreams.algo = 'twister';
            end
        end
        if ischar(a) && strcmpi(a,'reset')
            rand(DynareRandomStreams.algo,DynareRandomStreams.seed);
            randn('state',DynareRandomStreams.seed);
            return
        end
        if (~ischar(a) && isint(a)) || (ischar(a) && strcmpi(a,'clock'))
            if ischar(a)
                DynareRandomStreams.seed = floor(now*24*60*60);
            else
                DynareRandomStreams.seed = a;
            end
            rand(DynareRandomStreams.algo,DynareRandomStreams.seed);
            randn('state',DynareRandomStreams.seed);
            return
        end
        error('set_dynare_seed:: Something is wrong in the calling sequence!')
    else
        error('set_dynare_seed:: Cannot use more than one input argument with your version of MATLAB/Octave!')
    end
end
