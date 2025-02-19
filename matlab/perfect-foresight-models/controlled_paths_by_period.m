function controlled_paths_by_period = controlled_paths_by_period(M_, options_, info_period)
% Reorganize M_.perfect_foresight_controlled_paths period by period.
% “info_period” should only be set in a perfect foresight with expectation errors context,
% in which case only the information available at that period is extracted.

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

[periods, first_simulation_period] = get_simulation_periods(options_);

controlled_paths_by_period = struct('exogenize_id', cell(periods, 1), 'endogenize_id', cell(periods, 1), 'values', cell(periods, 1), 'learnt_in', cell(periods, 1));
for i=1:length(M_.perfect_foresight_controlled_paths)
    learnt_in = M_.perfect_foresight_controlled_paths(i).learnt_in;
    if isa(learnt_in, 'dates')
        learnt_in = learnt_in - first_simulation_period + 1;
    end
    if nargin >= 3 && learnt_in > info_period
        continue
    end

    prange = M_.perfect_foresight_controlled_paths(i).periods;
    if isa(prange, 'dates')
        prange = transpose(prange - first_simulation_period + 1);
    end
    if nargin >= 3 && any(prange < learnt_in)
        error('perfect_foresight_controlled_paths(learnt_in=%d): the periods are inconsistent with the learnt_in value', learnt_in)
    end
    for p = prange
        % Determine whether this is an change in expectations, and if yes overwrite the previous expected value
        idx = find(controlled_paths_by_period(p).exogenize_id == M_.perfect_foresight_controlled_paths(i).exogenize_id ...
                   & controlled_paths_by_period(p).endogenize_id == M_.perfect_foresight_controlled_paths(i).endogenize_id ...
                   & controlled_paths_by_period(p).learnt_in < learnt_in);
        if isempty(idx)
            controlled_paths_by_period(p).exogenize_id(end+1) = M_.perfect_foresight_controlled_paths(i).exogenize_id;
            controlled_paths_by_period(p).endogenize_id(end+1) = M_.perfect_foresight_controlled_paths(i).endogenize_id;
            controlled_paths_by_period(p).values(end+1) = M_.perfect_foresight_controlled_paths(i).value;
            controlled_paths_by_period(p).learnt_in(end+1) = learnt_in;
        else
            controlled_paths_by_period(p).values(idx) = M_.perfect_foresight_controlled_paths(i).value;
            controlled_paths_by_period(p).learnt_in(idx) = learnt_in;
        end
    end
end


%% Do various sanity checks

for p = 1:periods
    exogenize_id = controlled_paths_by_period(p).exogenize_id;
    [~, idx] = unique(exogenize_id, 'stable');
    duplicate_idx = setdiff(1:numel(exogenize_id), idx);
    if ~isempty(duplicate_idx)
        error('perfect_foresight_controlled_paths: variable %s is exogenized two times in period %d', M_.endo_names{exogenize_id(duplicate_idx(1))}, p);
    end

    endogenize_id = controlled_paths_by_period(p).endogenize_id;
    [~, idx] = unique(endogenize_id, 'stable');
    duplicate_idx = setdiff(1:numel(endogenize_id), idx);
    if ~isempty(duplicate_idx)
        error('perfect_foresight_controlled_paths: variable %s is endogenized two times in period %d', M_.exo_names{endogenize_id(duplicate_idx(1))}, p);
    end
end

if isfield(M_, 'det_shocks')
    for i = 1:length(M_.det_shocks)
        if M_.det_shocks(i).exo_det
            continue
        end
        prange = M_.det_shocks(i).periods;
        if isa(prange, 'dates')
            prange = transpose(prange - first_simulation_period + 1);
        end
        for p = prange
            exo_id = M_.det_shocks(i).exo_id;
            if ismember(exo_id, controlled_paths_by_period(p).endogenize_id)
                error('Exogenous variable %s is present in both the shocks and the perfect_foresight_controlled_paths in period %d', M_.exo_names{exo_id}, p);
            end
        end
    end
end

if nargin >= 3 && ~isempty(M_.learnt_shocks)
    for i = 1:length(M_.learnt_shocks)
        learnt_in = M_.learnt_shocks.learnt_in;
        if isa(learnt_in, 'dates')
            learnt_in = learnt_in - first_simulation_period + 1;
        end
        if learnt_in > info_period
            continue
        end

        prange = M_.learnt_shocks(i).periods;
        if isa(prange, 'dates')
            prange = transpose(prange - first_simulation_period + 1);
        end
        for p = prange
            exo_id = M_.learnt_shocks(i).exo_id;
            if ismember(exo_id, controlled_paths_by_period(p).endogenize_id)
                error('Exogenous variable %s is present in both the shocks(learnt_in=%d) and the perfect_foresight_controlled_paths in period %d', learnt_in, M_.exo_names{exo_id}, p);
            end
        end
    end
end
