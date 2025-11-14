function oo_ = correlation_mc_analysis(SampleSize,type,dname,fname,vartan,nvar,var1,var2,nar,mh_conf_sig,oo_,M_,options_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables correlation function.

% Copyright © 2008-2025 Dynare Team
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

if strcmpi(type,'posterior')
    folder_name=get_posterior_folder_name(options_);
    TYPE = 'Posterior';
    PATH = [dname filesep folder_name filesep];
else
    TYPE = 'Prior';
    PATH = [dname '/prior/moments/'];
end

ListOfFiles = dir([ PATH  fname '_' TYPE 'Correlations*.mat']);

[hh_fig, length_of_old_string] = wait_bar.run(0, [], 'Endogenous moments: correlation', options_.console_mode, 0);
if ~options_.console_mode
    set(hh_fig,'Name', 'Endogenous moments: correlation.' );
end

for var_iter_1=1:nvar
    i1 = 1; tmp = zeros(SampleSize,nvar,nvar,nar);
    for file = 1:length(ListOfFiles)
        load([ PATH  ListOfFiles(file).name ],'Correlation_array');
        i2 = i1 + rows(Correlation_array) - 1;
        tmp(i1:i2,:,:,:) = Correlation_array;
        i1 = i2+1;
    end

    var1=vartan{var_iter_1};
    for var_iter_2=1:nvar
        var2=vartan{var_iter_2};
        [hh_fig, length_of_old_string] = wait_bar.run(((var_iter_1-1)*nvar+var_iter_2)/(nvar^2), hh_fig, 'Endogenous moments: correlation', options_.console_mode, length_of_old_string);

        oo_ = initialize_output_structure(var1,var2,nar,TYPE,oo_,options_);
        for nar_iter=1:nar
            if options_.estimation.moments_posterior_density.indicator
                [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
                    posterior_moments(tmp(:,var_iter_1,var_iter_2,nar_iter),mh_conf_sig);
            else
                [p_mean, p_median, p_var, hpd_interval, p_deciles] = ...
                    posterior_moments(tmp(:,var_iter_1,var_iter_2,nar_iter),mh_conf_sig);
            end
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Mean',nar_iter,p_mean);
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Median',nar_iter,p_median);
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'Variance',nar_iter,p_var);
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDinf',nar_iter,hpd_interval(1));
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'HPDsup',nar_iter,hpd_interval(2));
            oo_ = fill_output_structure(var1,var2,TYPE,oo_,'deciles',nar_iter,p_deciles);
            if options_.estimation.moments_posterior_density.indicator
                oo_ = fill_output_structure(var1,var2,TYPE,oo_,'density',nar_iter,density);
            end
        end
    end
end
wait_bar.close(hh_fig,options_.console_mode)

function oo_ = initialize_output_structure(var1,var2,nar,type,oo_,options_)
oo_.([type, 'TheoreticalMoments']).dsge.correlation.Mean.(var1).(var2) = NaN(nar,1);
oo_.([type, 'TheoreticalMoments']).dsge.correlation.Median.(var1).(var2) = NaN(nar,1);
oo_.([type, 'TheoreticalMoments']).dsge.correlation.Variance.(var1).(var2) = NaN(nar,1);
oo_.([type, 'TheoreticalMoments']).dsge.correlation.HPDinf.(var1).(var2) = NaN(nar,1);
oo_.([type, 'TheoreticalMoments']).dsge.correlation.HPDsup.(var1).(var2) = NaN(nar,1);
oo_.([type, 'TheoreticalMoments']).dsge.correlation.deciles.(var1).(var2) = cell(nar,1);
if options_.estimation.moments_posterior_density.indicator
    oo_.([type, 'TheoreticalMoments']).dsge.correlation.density.(var1).(var2) = cell(nar,1);
end
for i=1:nar
    if options_.estimation.moments_posterior_density.indicator
        oo_.([type, 'TheoreticalMoments']).dsge.correlation.density.(var1).(var2)(i,1) = {NaN};
    end
    oo_.([type, 'TheoreticalMoments']).dsge.correlation.deciles.(var1).(var2)(i,1) = {NaN};
end

function oo_ = fill_output_structure(var1,var2,type,oo_,moment,lag,result)
switch moment
    case {'Mean','Median','Variance','HPDinf','HPDsup'}
        oo_.([type,  'TheoreticalMoments']).dsge.correlation.(moment).(var1).(var2)(lag,1) = result;
    case {'deciles','density'}
        oo_.([type, 'TheoreticalMoments']).dsge.correlation.(moment).(var1).(var2)(lag,1) = {result};
    otherwise
        disp('fill_output_structure:: Unknown field!')
end