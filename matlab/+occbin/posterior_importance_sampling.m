function oo_ = posterior_importance_sampling(M_, estim_params_, oo_, options_, bayestopt_)
% oo_ = posterior_importance_sampling(M_, estim_params_, oo_, options_, bayestopt_)
%
% Uses draws from previous MCMC (either linear of PKF) 
% to re-draw parameters using PKF/PPF with importance sampling
%
% INPUTS
%  - M_                     [structure] MATLAB's structure describing the model
%  - estim_params_          [structure] structure containing estimation parameters information
%  - oo_                    [structure] MATLAB's structure storing the results
%  - options_               [structure] MATLAB's structure containing the options
%  - bayestopt_             [structure] MATLAB's structure containing Bayesian estimation options
%
% OUTPUTS
%  - oo_                    [structure] MATLAB's structure storing the results with OccBin posterior moments and statistics

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

DirectoryName = CheckPath('metropolis',M_.dname);
OutputDirectoryName = CheckPath('Output',M_.dname);
% retrieve MCMC from original estimation
orig_dname = options_.occbin.posterior_importance_sampling.orig_dname;
orig_fname = options_.occbin.posterior_importance_sampling.orig_fname;
OrigDirectoryName = CheckPath('metropolis',orig_dname);
a=dir([OrigDirectoryName '/'  orig_fname '_mh_history*']);
load([OrigDirectoryName '/'  orig_fname '_mh_history_' int2str(length(a)-1)]);
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
npar = length(bayestopt_.name);

KeptNumberOfMhDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
 
B = options_.sub_draws;

if B==KeptNumberOfMhDraws*options_.mh_nblck
    % we load all retained MH runs !
    lpost = GetAllPosteriorDraws(options_, orig_dname, orig_fname, 0, record.KeepedDraws.FirstLine, record.KeepedDraws.FirstMhFile);
    params0 = GetAllPosteriorDraws(options_, orig_dname, orig_fname, 'all', record.KeepedDraws.FirstLine, record.KeepedDraws.FirstMhFile);

else
    M_bkp = M_;
    M_.dname = orig_dname;
    M_.fname = orig_fname;
    [params0, lpost]=get_posterior_subsample(M_,options_,B);    
    M_ = M_bkp;
end

% compute posterior kernel with occbin and derive weights
lpk = nan(size(lpost));
options_.occbin.likelihood.status=true;
this_filter = options_.occbin.posterior_importance_sampling.filter;
orig_filter = options_.occbin.posterior_importance_sampling.orig_filter;
if isequal(this_filter, 'ppf')
    options_.occbin.filter.particle.status=true;
elseif isequal(this_filter, 'pkf')
    options_.occbin.filter.particle.status=false;
    options_.occbin.likelihood.IF_likelihood=false; % make sure no IF
else
    error('invalid filter name')
end
% make sure nothing is plotted during sampling
options_.occbin.filter.particle.nograph=true; 
options_.occbin.filter.particle.likelihood_only=true;

Label = '%s %s';
bar_title = 'Posterior importance sampling';
[hh_fig, length_of_old_string] = wait_bar.run(0, [], sprintf(Label, bar_title, '...'), options_.console_mode, 0, bar_title, 0);
bar_title = ['Posterior kernel with ' this_filter ' filter'];
waitbar_title = sprintf(Label, bar_title, ' %s');

for k=1:length(lpk)
    lpk(k) = evaluate_posterior_kernel(params0(k,:)', M_, estim_params_, oo_, options_, bayestopt_);
    txt = sprintf('Step %u of %u', k,length(lpk));
    [~, length_of_old_string]=wait_bar.run(k/length(lpk), hh_fig, sprintf(waitbar_title, txt),options_.console_mode,length_of_old_string,[],0);
    
end
wait_bar.close(hh_fig,options_.console_mode);

lnw0 = lpk-lpost;
lnw0 = lnw0-max(lnw0);
weights0 = ones(1,length(lpk))/length(lpk) ;
wtilde0 = weights0.*exp(lnw0);
weights = wtilde0/sum(wtilde0);
ESS = 1/sum(weights.^2);

oo_.occbin.ESS = ESS;

fprintf(' The effective sample size (ESS) is %.2f for %.0f draws\n',ESS,B)

% weighted posterior bootstrap sample with OccBin posterior kernel
resampled_params = resample(params0,weights,options_.particle);
save([DirectoryName '/'  M_.fname '_occbin_importance_sampling'],'params0','resampled_params','lpost','lpk','ESS','weights')

% get mh mode of OccBin posterior kernel
[~, im] = max(lpk);
xparam1 = params0(im,:)';
fval = lpk(im);
hh = inv(cov(resampled_params));
parameter_names = bayestopt_.name;
save([OutputDirectoryName '/'  M_.fname '_mhocc_mode'],'fval','hh','parameter_names','xparam1')

% print results and store them
nvx     = estim_params_.nvx;
nvn     = estim_params_.nvn;
ncx     = estim_params_.ncx;
ncn     = estim_params_.ncn;
np      = estim_params_.np ;
nxs     = nvx+nvn+ncx+ncn;

oo1 = oo_;

skipline(2)
disp('ESTIMATION RESULTS: OccBin importance sampling from linear MCMC')
skipline()

header_width = row_header_width(M_, estim_params_, bayestopt_);

hpd_interval=[num2str(options_.mh_conf_sig*100), '% HPD interval'];
tit2 = sprintf('%-*s %8s %12s %12s %23s %8s %12s\n',header_width,' ','filter','prior mean','post. mean',hpd_interval,'prior','pstdev');
pformat =      '%-*s %-8s %12.3f %12.4f %11.4f %11.4f %8s %12.4f';
pformat1 =     '%-*s %-20s  %12.4f %11.4f %11.4f';
pnames = prior_dist_names;

if np
    type = 'parameters';
    skipline()
    disp(type)
    disp(tit2)

    for i=1:np
        ip = i+nxs;
        Draws = resampled_params(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        name = bayestopt_.name{ip};
        oo1 = Filloo(oo1, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat, header_width, name, this_filter, bayestopt_.p1(ip),...
            post_mean, ...
            hpd_interval, ...
            pnames{bayestopt_.pshape(ip)+1}, ...
            bayestopt_.p2(ip)));
        Draws = params0(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat1, header_width, ' ', orig_filter,...
            post_mean, ...
            hpd_interval));
    end
end

if nvx
    type = 'shocks_std';
    skipline()
    disp('standard deviation of shocks')
    disp(tit2)
    for ip=1:nvx
        Draws = resampled_params(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        name = bayestopt_.name{ip};
        oo1 = Filloo(oo1, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat, header_width, name, this_filter, bayestopt_.p1(ip),...
            post_mean, ...
            hpd_interval, ...
            pnames{bayestopt_.pshape(ip)+1}, ...
            bayestopt_.p2(ip)));
        Draws = params0(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat1, header_width, ' ', orig_filter,...
            post_mean, ...
            hpd_interval));
    end
end

if nvn
    type = 'measurement_errors_std';
    skipline()
    disp('standard deviation of measurement errors')
    disp(tit2)
    for i=1:nvn
        ip=i+nvx;
        Draws = resampled_params(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        name = bayestopt_.name{ip};
        oo1 = Filloo(oo1, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat, header_width, name, this_filter, bayestopt_.p1(ip),...
            post_mean, ...
            hpd_interval, ...
            pnames{bayestopt_.pshape(ip)+1}, ...
            bayestopt_.p2(ip)));
        Draws = params0(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat1, header_width, ' ', orig_filter,...
            post_mean, ...
            hpd_interval));
    end
end

if ncx
    type = 'shocks_corr';
    skipline()
    disp('correlation of shocks')
    disp(tit2)
    for i=1:ncx
        ip=i+nvx+nvn;
        Draws = resampled_params(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        name = bayestopt_.name{ip};
        name = strrep(name,'corr ','');
        NAME = strrep(name,', ','_');
        name = strrep(name,', ',',');
        oo1 = Filloo(oo1, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat, header_width, name, this_filter, bayestopt_.p1(ip),...
            post_mean, ...
            hpd_interval, ...
            pnames{bayestopt_.pshape(ip)+1}, ...
            bayestopt_.p2(ip)));
        Draws = params0(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat1, header_width, ' ', orig_filter,...
            post_mean, ...
            hpd_interval));
    end
end

if ncn
    type = 'measurement_errors_corr';
    skipline()
    disp('correlation of measurement errors')
    disp(tit2)
    for i=1:ncn
        ip=i+nvx+nvn+ncx;
        Draws = resampled_params(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        name = bayestopt_.name{ip};
        name = strrep(name,'corr ','');
        NAME = strrep(name,', ','_');
        name = strrep(name,', ',',');
        oo1 = Filloo(oo1, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat, header_width, name, this_filter, bayestopt_.p1(ip),...
            post_mean, ...
            hpd_interval, ...
            pnames{bayestopt_.pshape(ip)+1}, ...
            bayestopt_.p2(ip)));
        Draws = params0(:,ip);
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws, options_.mh_conf_sig);
        oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        disp(sprintf(pformat1, header_width, ' ', orig_filter,...
            post_mean, ...
            hpd_interval));
    end
end

% this  is where we store posterior statistics after importance sampling
field_list = {'posterior_mean', 'posterior_hpdinf', 'posterior_hpdsup', 'posterior_median', 'posterior_variance', 'posterior_std', 'posterior_deciles', 'posterior_density'};
for k=1:length(field_list)
    oo_.occbin.(field_list{k}) = oo1.(field_list{k});
end

if not(options_.nograph)
    OutputDirectoryName = CheckPath('graphs',M_.dname);
    TeX     = false; %options_.TeX;
    MaxNumberOfPlotPerFigure = 9;% The square root must be an integer!
    nn = sqrt(MaxNumberOfPlotPerFigure);

    figurename = [this_filter ' and ' orig_filter ' posteriors'];
    figunumber = 0;
    subplotnum = 0;
    for i=1:npar
        subplotnum = subplotnum+1;
        if subplotnum == 1
            figunumber = figunumber+1;
            hfig=dyn_figure(options_.nodisplay, 'Name', figurename);
        end
        [nam,texnam] = get_the_name(i, TeX, M_, estim_params_, options_);
        if i <= nvx
            name = M_.exo_names{estim_params_.var_exo(i,1)};
            x1 = oo_.occbin.posterior_density.shocks_std.(name)(:,1);
            f1 = oo_.occbin.posterior_density.shocks_std.(name)(:,2);
            x2 = oo_.posterior_density.shocks_std.(name)(:,1);
            f2 = oo_.posterior_density.shocks_std.(name)(:,2);
            if ~options_.mh_posterior_mode_estimation
                pmod = oo_.posterior_mode.shocks_std.(name);
            end
        elseif i <= nvx+nvn
            name = options_.varobs{estim_params_.nvn_observable_correspondence(i-nvx,1)};
            x1 = oo_.occbin.posterior_density.measurement_errors_std.(name)(:,1);
            f1 = oo_.occbin.posterior_density.measurement_errors_std.(name)(:,2);
            x2 = oo_.posterior_density.measurement_errors_std.(name)(:,1);
            f2 = oo_.posterior_density.measurement_errors_std.(name)(:,2);
            if ~options_.mh_posterior_mode_estimation
                pmod = oo_.posterior_mode.measurement_errors_std.(name);
            end
        elseif i <= nvx+nvn+ncx
            j = i - (nvx+nvn);
            k1 = estim_params_.corrx(j,1);
            k2 = estim_params_.corrx(j,2);
            name = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
            x1 = oo_.occbin.posterior_density.shocks_corr.(name)(:,1);
            f1 = oo_.occbin.posterior_density.shocks_corr.(name)(:,2);
            x2 = oo_.posterior_density.shocks_corr.(name)(:,1);
            f2 = oo_.posterior_density.shocks_corr.(name)(:,2);
            if ~options_.mh_posterior_mode_estimation
                pmod = oo_.posterior_mode.shocks_corr.(name);
            end
        elseif i <= nvx+nvn+ncx+ncn
            j = i - (nvx+nvn+ncx);
            k1 = estim_params_.corrn(j,1);
            k2 = estim_params_.corrn(j,2);
            name = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
            x1 = oo_.occbin.posterior_density.measurement_errors_corr.(name)(:,1);
            f1 = oo_.occbin.posterior_density.measurement_errors_corr.(name)(:,2);
            x2 = oo_.posterior_density.measurement_errors_corr.(name)(:,1);
            f2 = oo_.posterior_density.measurement_errors_corr.(name)(:,2);
            if ~options_.mh_posterior_mode_estimation
                pmod = oo_.posterior_mode.measurement_errors_corr.(name);
            end
        else
            j = i - (nvx+nvn+ncx+ncn);
            name = M_.param_names{estim_params_.param_vals(j,1)};
            x1 = oo_.occbin.posterior_density.parameters.(name)(:,1);
            f1 = oo_.occbin.posterior_density.parameters.(name)(:,2);
            x2 = oo_.posterior_density.parameters.(name)(:,1);
            f2 = oo_.posterior_density.parameters.(name)(:,2);
            if ~options_.mh_posterior_mode_estimation
                pmod = oo_.posterior_mode.parameters.(name);
            end
        end
        top1 = max(f1);
        top2 = max(f2);
        top0 = max([top1; top2]);
        binf1 = x1(1);
        bsup1 = x1(end);
        binf2 = x2(1);
        bsup2 = x2(end);
        borneinf = min(binf1, binf2);
        bornesup = max(bsup1, bsup2);
        subplot(nn, nn, subplotnum)
        plot(x2, f2, '-b', 'linewidth', 2);
        hold on;
        plot(x1, f1, '--r', 'linewidth', 2);
        box on
        axis([borneinf bornesup 0 1.1*top0])
        if TeX
            title(texnam, 'Interpreter', 'latex')
        else
            title(nam, 'Interpreter', 'none')
        end
        hold off
        drawnow
        if subplotnum == MaxNumberOfPlotPerFigure || i == npar
            dyn_saveas(hfig,[OutputDirectoryName '/' M_.fname '_PosteriorsOccbinIS_' int2str(figunumber)], options_.nodisplay, options_.graph_format);
            if TeX && any(strcmp('eps', cellstr(options_.graph_format)))
                fprintf(fidTeX, '\\begin{figure}[H]\n');
                fprintf(fidTeX, '\\centering\n');
                fprintf(fidTeX, '\\includegraphics[width=%2.2f\\textwidth]{%s/%s_PosteriorsOccbinIS_%s}\n', ...
                    options_.figures.textwidth*min(subplotnum/nn,1), OutputDirectoryName, M_.fname, int2str(figunumber));
                fprintf(fidTeX,['\\caption{' this_filter ' and ' orig_filter ' posteriors.}']);
                fprintf(fidTeX,'\\label{Fig:PosteriorsOccbinIS:%s}\n', int2str(figunumber));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
                if i == npar
                    fprintf(fidTeX,'%% End of TeX file.\n');
                    fclose(fidTeX);
                end
            end
            subplotnum = 0;
        end
    end
end

end

function oo = Filloo(oo, name, type, postmean, hpdinterval, postmedian, postvar, postdecile, density)
oo.posterior_mean.(type).(name) = postmean;
oo.posterior_hpdinf.(type).(name) = hpdinterval(1);
oo.posterior_hpdsup.(type).(name) = hpdinterval(2);
oo.posterior_median.(type).(name) = postmedian;
oo.posterior_variance.(type).(name) = postvar;
oo.posterior_std.(type).(name) = sqrt(postvar);
oo.posterior_deciles.(type).(name) = postdecile;
oo.posterior_density.(type).(name) = density;

end