function oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,table_title,field_name)
%function oo_=display_results_table(xparam1,stdh,M_,estim_params_,bayestopt_,oo_,pnames,table_title,field_name)
% Display estimation results on screen and write them to TeX-file
%
% INPUTS
%   o xparam1       [double]   (p*1) vector of estimate parameters.
%   o stdh          [double]   (p*1) vector of estimate parameters.
%   o M_                        MATLAB's structure describing the Model
%   o estim_params_             MATLAB's structure describing the estimated_parameters
%   o options_                  MATLAB's structure describing the options
%   o bayestopt_                MATLAB's structure describing the priors
%   o oo_                       MATLAB's structure gathering the results
%   o pnames        [string]    Cell of strings storing the names for prior distributions
%   o table_title   [string]    Title of the Table
%   o field_name    [string]    String storing the name of the fields for oo_ where the parameters are stored
%
% OUTPUTS
%   o oo_                       MATLAB's structure gathering the results
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2014-2025 Dynare Team
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

skipline()
disp(['RESULTS FROM ' upper(table_title) ' ESTIMATION'])
LaTeXtitle=strrep(table_title,' ','_');
tstath = abs(xparam1)./stdh;

header_width = row_header_width(M_, estim_params_, bayestopt_);
if contains(field_name,'posterior')
    tit1 = sprintf('%-*s %10s %8s %7s %6s %6s\n', header_width, ' ', 'prior mean', ...
                    'mode', 's.d.', 'prior', 'pstdev');
else
    tit1 = sprintf('%-*s %10s %7s %6s\n', header_width, ' ', 'Estimate', 's.d.', 't-stat');
end

if estim_params_.np % number of estimated structural parameters
    ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+1; % offset: structural parameters are ordered last in xparam1
    disp('parameters')
    disp(tit1)
    for i=1:estim_params_.np
        name = bayestopt_.name{ip};
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width,name, ...
                    bayestopt_.p1(ip),xparam1(ip),stdh(ip), ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n', ...
                    header_width, name, xparam1(ip), stdh(ip), tstath(ip));
        end
        oo_.(sprintf('%s_mode', field_name)).parameters.(name) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).parameters.(name) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if estim_params_.nvx % number of estimated stderr parameters for structural shocks
    ip = 1; % offset: stderr of shocks are ordered first in xparam1
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:estim_params_.nvx
        k = estim_params_.var_exo(i,1);
        name = M_.exo_names{k};
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width, name, bayestopt_.p1(ip), xparam1(ip), ...
                    stdh(ip), pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n', header_width, name, xparam1(ip), stdh(ip), tstath(ip));
        end
        M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
        oo_.(sprintf('%s_mode', field_name)).shocks_std.(name) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).shocks_std.(name) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if estim_params_.nvn % number of estimated stderr parameters for measurement errors
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = estim_params_.nvx+1; % offset: stderr of measurement errors are ordered second in xparam1
    for i=1:estim_params_.nvn
        name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width, name, bayestopt_.p1(ip), ...
                    xparam1(ip), stdh(ip), ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n', header_width, name, xparam1(ip), ...
                    stdh(ip), tstath(ip))
        end
        oo_.(sprintf('%s_mode', field_name)).measurement_errors_std.(name) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).measurement_errors_std.(name) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if estim_params_.ncx % number of estimated corr parameters for structural shocks
    disp('correlation of shocks')
    disp(tit1)
    ip = estim_params_.nvx+estim_params_.nvn+1; % offset: corr of shocks are ordered third in xparam1
    for i=1:estim_params_.ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        name = sprintf('%s,%s', M_.exo_names{k1}, M_.exo_names{k2});
        NAME = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width, name, bayestopt_.p1(ip), xparam1(ip), stdh(ip),  ...
                    pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n', header_width,name, xparam1(ip), ...
                    stdh(ip), tstath(ip));
        end
        M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
        M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
        oo_.(sprintf('%s_mode', field_name)).shocks_corr.(NAME) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).shocks_corr.(NAME) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if estim_params_.ncn % number of estimated corr parameters for measurement errors
    disp('correlation of measurement errors')
    disp(tit1)
    ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+1; % offset: corr of measurement errors are ordered fourth in xparam1
    for i=1:estim_params_.ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        name = sprintf('%s,%s', M_.endo_names{k1}, M_.endo_names{k2});
        NAME = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width, name, bayestopt_.p1(ip), xparam1(ip), stdh(ip), ...
                    pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n',header_width, name, xparam1(ip), ...
                    stdh(ip), tstath(ip));
        end
        oo_.(sprintf('%s_mode', field_name)).measurement_errors_corr.(NAME) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).measurement_errors_corr.(NAME) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if estim_params_.nsx % number of estimated skew parameters for structural shocks
    disp('skewness of shocks')
    disp(tit1)
    ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+1; % offset: skew of shocks are ordered fifth in xparam1
    for i=1:estim_params_.nsx
        k = estim_params_.skew_exo(i,1);
        name = sprintf('%s', M_.exo_names{k});
        NAME = sprintf('%s', M_.exo_names{k});
        if contains(field_name,'posterior')
            fprintf('%-*s %10.4f %8.4f %7.4f %6s %6.4f \n', ...
                    header_width, name, bayestopt_.p1(ip), xparam1(ip), stdh(ip), ...
                    pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        else
            fprintf('%-*s %10.4f %7.4f %7.4f \n',header_width, name, xparam1(ip), ...
                    stdh(ip), tstath(ip));
        end
        oo_.(sprintf('%s_mode', field_name)).shocks_skew.(NAME) = xparam1(ip);
        oo_.(sprintf('%s_std_at_mode', field_name)).shocks_skew.(NAME) = stdh(ip);
        ip = ip+1;
    end
    skipline()
end

if any(xparam1(1:estim_params_.nvx+estim_params_.nvn)<0)
    warning(sprintf('Some estimated standard deviations are negative.\n         Dynare internally works with variances so that the sign does not matter.\n         Nevertheless, it is recommended to impose either prior restrictions (Bayesian Estimation)\n         or a lower bound (ML) to assure positive values.'))
end

latexDirectoryName = CheckPath('latex',M_.dname);

if any(bayestopt_.pshape > 0) && options_.TeX %% Bayesian estimation (posterior mode) Latex output
    if estim_params_.np % structural parameters
        filename = [latexDirectoryName '/' M_.fname '_Posterior_Mode_1.tex'];
        fidTeX = fopen(filename,'w');
        TeXBegin_Bayesian(fidTeX,1,'parameters')
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+1; % offset: structural parameters are ordered last in xparam1
        for i=1:estim_params_.np
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    M_.param_names_tex{estim_params_.param_vals(i,1)}, ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip + 1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nvx % stderr of shocks
        TeXfile = [latexDirectoryName '/' M_.fname '_Posterior_Mode_2.tex'];
        fidTeX = fopen(TeXfile,'w');
        TeXBegin_Bayesian(fidTeX,2,'standard deviation of structural shocks')
        ip = 1; % offset: stderr of shocks are ordered first in xparam1
        for i=1:estim_params_.nvx
            k = estim_params_.var_exo(i,1);
            fprintf(fidTeX, '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    M_.exo_names_tex{k},...
                    pnames{bayestopt_.pshape(ip)+1},...
                    bayestopt_.p1(ip),...
                    bayestopt_.p2(ip),...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nvn % stderr of measurement errors
        TeXfile = [latexDirectoryName '/' M_.fname '_Posterior_Mode_3.tex'];
        fidTeX  = fopen(TeXfile,'w');
        TeXBegin_Bayesian(fidTeX,3,'standard deviation of measurement errors')
        ip = estim_params_.nvx+1; % offset: stderr of measurement errors are ordered second in xparam1
        for i=1:estim_params_.nvn
            idx = strmatch(options_.varobs{estim_params_.nvn_observable_correspondence(i,1)}, M_.endo_names);
            fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    M_.endo_names_tex{idx}, ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip),...
                    xparam1(ip),...
                    stdh(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.ncx % corr of shocks
        TeXfile = [latexDirectoryName '/' M_.fname '_Posterior_Mode_4.tex'];
        fidTeX = fopen(TeXfile,'w');
        TeXBegin_Bayesian(fidTeX,4,'correlation of structural shocks')
        ip = estim_params_.nvx+estim_params_.nvn+1; % offset: corr of shocks are ordered third in xparam1
        for i=1:estim_params_.ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            fprintf(fidTeX, '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    [M_.exo_names_tex{k1} ',' M_.exo_names_tex{k2}], ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.ncn % corr of measurement errors
        TeXfile = [latexDirectoryName '/' M_.fname '_Posterior_Mode_5.tex'];
        fidTeX = fopen(TeXfile,'w');
        TeXBegin_Bayesian(fidTeX,5,'correlation of measurement errors')
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+1; % offset: corr of measurement errors are ordered fourth in xparam1
        for i=1:estim_params_.ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    [ M_.endo_names_tex{k1} ',' M_.endo_names_tex{k2}], ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nsx % skew of shocks
        TeXfile = [latexDirectoryName '/' M_.fname '_Posterior_Mode_6.tex'];
        fidTeX = fopen(TeXfile,'w');
        TeXBegin_Bayesian(fidTeX,6,'skeweness of structural shocks')
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+1; % offset: skew of shocks are ordered fifth in xparam1
        for i=1:estim_params_.nsx
            k = estim_params_.skew_exo(i,1);
            name = [M_.exo_names_tex{k}];
            fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    name, ...
                    pnames{bayestopt_.pshape(ip)+1}, ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
elseif all(bayestopt_.pshape == 0) && options_.TeX %% MLE and GMM Latex output
    if estim_params_.np % structural parameters
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_1.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 1, 'parameters', table_title, LaTeXtitle)
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+1; % offset: structural parameters are ordered last in xparam1
        for i=1:estim_params_.np
            fprintf(fidTeX,'$%s$ & %8.4f & %7.4f & %7.4f\\\\ \n',...
                    M_.param_names_tex{estim_params_.param_vals(i,1)}, ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip + 1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nvx % stderr of shocks
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_2.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 2, 'standard deviation of structural shocks', table_title, LaTeXtitle)
        ip = 1; % offset: stderr of shocks are ordered first in xparam1
        for i=1:estim_params_.nvx
            k = estim_params_.var_exo(i,1);
            fprintf(fidTeX, '$%s$ & %8.4f & %7.4f & %7.4f\\\\ \n', ...
                    M_.exo_names_tex{k}, ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nvn % stderr of measurement errors
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_3.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 3, 'standard deviation of measurement errors', table_title, LaTeXtitle)
        ip = estim_params_.nvx+1; % offset: stderr of measurement errors are ordered second in xparam1
        for i=1:estim_params_.nvn
            idx = strmatch(options_.varobs{estim_params_.nvn_observable_correspondence(i,1)}, M_.endo_names);
            fprintf(fidTeX, '$%s$ & %8.4f & %7.4f & %7.4f \\\\ \n', ...
                    M_.endo_names_tex{idx}, ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.ncx % corr of shocks
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_4.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 4, 'correlation of structural shocks', table_title,LaTeXtitle)
        ip = estim_params_.nvx+estim_params_.nvn+1; % offset: corr of shocks are ordered third in xparam1
        for i=1:estim_params_.ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            fprintf(fidTeX, '$%s$  & %8.4f & %7.4f & %7.4f \\\\ \n', ...
                    [M_.exo_names_tex{k1} ',' M_.exo_names_tex{k2}], ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.ncn % corr of measurement errors
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_5.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 5, 'correlation of measurement errors', table_title, LaTeXtitle)
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+1; % offset: corr of measurement errors are ordered fourth in xparam1
        for i=1:estim_params_.ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            fprintf(fidTeX, '$%s$  & %8.4f & %7.4f & %7.4f \\\\ \n', ...
                    [ M_.endo_names_tex{k1} ',' M_.endo_names_tex{k2}], ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
    if estim_params_.nsx % skew of shocks
        filename = [latexDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_6.tex'];
        fidTeX = fopen(filename, 'w');
        TeXBegin_ML(fidTeX, 6, 'skewness of structural shocks', table_title, LaTeXtitle)
        ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+1; % offset: skew of shocks are ordered fifth in xparam1
        for i=1:estim_params_.nsx
            k = estim_params_.skew_exo(i,1);
            name = [M_.exo_names_tex{k}];
            fprintf(fidTeX, '$%s$  & %8.4f & %7.4f & %7.4f \\\\ \n', ...
                    name, ...
                    xparam1(ip), ...
                    stdh(ip), ...
                    tstath(ip));
            ip = ip+1;
        end
        TeXEnd(fidTeX)
    end
end



%% subfunctions:
%
function TeXBegin_Bayesian(fid, fnum, title)
fprintf(fid,'%% TeX-table generated by dynare_estimation (Dynare).\n');
fprintf(fid,['%% RESULTS FROM POSTERIOR MAXIMIZATION (' title ')\n']);
fprintf(fid,['%% ' datestr(now,0)]);
fprintf(fid,' \n');
fprintf(fid,' \n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{longtable}{llcccc} \n');
fprintf(fid,['\\caption{Results from posterior maximization (' title ')}\\\\\n ']);
fprintf(fid,['\\label{Table:Posterior:' int2str(fnum)  '}\\\\\n']);
fprintf(fid,'\\toprule \n');
fprintf(fid,'  & \\multicolumn{3}{c}{Prior}  &  \\multicolumn{2}{c}{Posterior} \\\\\n');
fprintf(fid,'  \\cmidrule(r{.75em}){2-4} \\cmidrule(r{.75em}){5-6}\n');
fprintf(fid,'  & Dist. & Mean  & Stdev & Mode & Stdev \\\\ \n');
fprintf(fid,'\\midrule \\endfirsthead \n');
fprintf(fid,'\\caption{(continued)}\\\\\n ');
fprintf(fid,'\\bottomrule \n');
fprintf(fid,'  & \\multicolumn{3}{c}{Prior}  &  \\multicolumn{2}{c}{Posterior} \\\\\n');
fprintf(fid,'  \\cmidrule(r{.75em}){2-4} \\cmidrule(r{.75em}){5-6}\n');
fprintf(fid,'  & Dist. & Mean  & Stdev & Mode & Stdev \\\\ \n');
fprintf(fid,'\\midrule \\endhead \n');
fprintf(fid,'\\bottomrule \\multicolumn{6}{r}{(Continued on next page)}\\endfoot \n');
fprintf(fid,'\\bottomrule\\endlastfoot \n');

function TeXBegin_ML(fid, fnum, title, table_title, LaTeXtitle)
fprintf(fid,'%% TeX-table generated by dynare_estimation (Dynare).\n');
fprintf(fid,['%% RESULTS FROM ' table_title ' MAXIMIZATION (' title ')\n']);
fprintf(fid,['%% ' datestr(now,0)]);
fprintf(fid,' \n');
fprintf(fid,' \n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{longtable}{llcc} \n');
fprintf(fid,['\\caption{Results from ' table_title ' maximization (' title ')}\\\\\n ']);
fprintf(fid,['\\label{Table:' LaTeXtitle ':' int2str(fnum) '}\\\\\n']);
fprintf(fid,'\\toprule \n');
fprintf(fid,'  & Mode & s.d. & t-stat\\\\ \n');
fprintf(fid,'\\midrule \\endfirsthead \n');
fprintf(fid,'\\caption{(continued)}\\\\\n ');
fprintf(fid,'\\toprule \n');
fprintf(fid,'  & Mode & s.d. & t-stat\\\\ \n');
fprintf(fid,'\\midrule \\endhead \n');
fprintf(fid,'\\bottomrule  \\multicolumn{4}{r}{(Continued on next page)} \\endfoot \n');
fprintf(fid,'\\bottomrule \\endlastfoot \n');

function TeXEnd(fid)
fprintf(fid,'\\end{longtable}\n ');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'%% End of TeX file.\n');
fclose(fid);
