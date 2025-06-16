function lprobs_sample = trace_plot_dime(options_, M_)

% Plot the history of the densities of an ensemble to visually inspect convergence.
%
% INPUTS
% - options_         [struct]   dynare's options
% - M_               [struct]   model description
% - oo_              [struct]   outputs
%
% SPECIAL REQUIREMENTS
% None.

% Copyright © 2022-2025 Dynare Team
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

lprobs = GetAllPosteriorDraws(options_, M_.dname, [], 0);
tune = options_.posterior_sampler_options.current_options.tune;
[niter, nchain] = size(lprobs);

if max(lprobs(end,:)) > 0
    ulim = max(lprobs(end,:))*1.2;
else
    ulim = max(lprobs(end,:))/1.2;
end
if min(lprobs(end,:)) > 0
    llim = min(lprobs(end,:))/5;
else
    llim = min(lprobs(end,:))*5;
end

graphFolder = CheckPath('graphs',M_.dname);
hh_fig = dyn_figure(options_.nodisplay,'Name','DIME Convergence Diagnostics');

hold on
lines_tune = plot(1:niter-tune, lprobs(1:end-tune,:), 'Color', '#D95319');
lines_sample = plot(niter-tune:niter, lprobs(end-tune:end,:), 'Color', '#0072BD');
if ~isoctave
    % Set the transparency (alpha channel) of line objects (undocumented MATLAB feature)
    % See https://fr.mathworks.com/matlabcentral/discussions/ideas/833472-add-alpha-capability-to-line-class
    for i = 1:nchain
        lines_tune(i).Color(4) = min(1,10/nchain);
        lines_sample(i).Color(4) = min(1,10/nchain);
    end
end
ylim([llim, ulim])
hold off

dyn_saveas(hh_fig,[graphFolder '/' M_.fname '_trace_lprob'],options_.nodisplay,options_.graph_format);
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_trace_lprob}\n',[graphFolder '/' M_.fname]);
    fprintf(fidTeX,'\\caption{Ensemble traces of posterior densities for DIME.\n');
    fprintf(fidTeX,'\\label{Fig:DIME_trace}\n');
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
end

lprobs_sample = lprobs(end-tune:end,:);
