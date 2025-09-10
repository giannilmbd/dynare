function [fig_handle, length_of_old_string]= run(prctdone, wait_handle, running_text, console_mode, move_cursor_backwards, fig_title, whoiam, Parallel)
% [fig_handle, length_of_old_string]= run(prctdone, wait_handle, running_text, console_mode, move_cursor_backwards, whoiam, Parallel)
% adaptive waitbar, producing console mode waitbars with
% GNU Octave and when console_mode=1
% Inputs
%  - prctdone               [double]    progress of waitbar
%  - wait_handle            [handle]    figure handle to waitbar (empty if parallel or console); not required for initialization, i.e. prctdone==0
%  - running_text           [char]      text displayed above waitbar
%  - console_mode           [bool]      indicator whether console output is requested
%  - move_cursor_backwards  [integer]   number of backward cursor movements for console
%  - fig_title              [string]    option figure title for initialization (prctdone==0)
%  - whoiam                 [integer]   In concurrent programming a modality to refer to the different threads running in parallel is needed.
%                                       The integer whoaim is the integer that allows us to distinguish between them. Then it is the 
%                                       index number of this CPU among all CPUs in the cluster.
%  - ThisMatlab             [integer]   Allows us to distinguish between the 'main' MATLAB, the slave MATLAB worker, local MATLAB, remote MATLAB,
%                                       then it is the index number of this slave machine in the cluster.
%
% Outputs
%  - fig_handle             [handle]    figure handle to waitbar (empty if parallel or console)
%  - length_of_old_string   [integer]   length of last text printed in console
%
% Copyright © 2011-2025 Dynare Team
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

if prctdone==0
    init=1;
else
    init=0;
end

if nargin< 6
    fig_title=[];
end
if nargin< 7 || isempty(whoiam) 
    whoiam=0;
end

% initialize outputs
fig_handle=[];
length_of_old_string=0;

if ~whoiam %serial execution
    if console_mode
        if init
            diary off;
            return
        end        
        newString=sprintf([running_text,' %3.f%% done'], prctdone*100);
        move_backwards=repmat('\b',1,move_cursor_backwards);
        fprintf([move_backwards,'%s'],newString);
        length_of_old_string=length(newString);
    else
        if init
            fig_handle = waitbar(prctdone,running_text);
            if ~isempty(fig_title)
                set(fig_handle,'Name', fig_title);
            end
        else
            fig_handle = waitbar(prctdone,wait_handle,running_text);
        end
    end
else %parallel execution
    if Parallel.Local
        waitbarTitle='Local ';
    else
        waitbarTitle=[Parallel.ComputerName];
    end
    fMessageStatus(prctdone,whoiam,running_text, waitbarTitle, Parallel);
end