function online_auxiliary_filter(xparam1, dataset_, options_, M_, estim_params_, bayestopt_, oo_)
% online_auxiliary_filter(xparam1, dataset_, options_, M_, estim_params_, bayestopt_, oo_)
% Liu & West particle filter = auxiliary particle filter including Liu & West filter on parameters.
%
% INPUTS
% - xparam1                  [double]    n×1 vector, Initial condition for the estimated parameters.
% - dataset_                 [dseries]   Sample used for estimation.
% - options_                 [struct]    Option values.
% - M_                       [struct]    Description of the model.
% - estim_params_            [struct]    Description of the estimated parameters.
% - bayestopt_               [struct]    Prior definition.
% - oo_                      [struct]    Results.
%
% Reference: Liu, Jane and Mike West (2001): “Combined parameter and state estimation in simulation-based filtering”,
% in Sequential Monte Carlo Methods in Practice, Eds. Doucet, Freitas and Gordon, Springer Verlag, Chapter
% 10, 197-223.

% Copyright © 2013-2026 Dynare Team
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

% Set seed for randn().
set_dynare_seed_local_options([],false,'default');
options_.verbosity=0; %particularly suppress warning messages during k_order_pert within the loop
pruning = options_.particle.pruning;
variance_update = true;
online_opt = options_.posterior_sampler_options.current_options;

% Set location for the simulated particles.
SimulationFolder = CheckPath('online', M_.dname);

bounds = prior_bounds(bayestopt_, options_.prior_trunc); % Reset bounds as lb and ub must only be operational during mode-finding

% initialization of state particles
[~, ~, ReducedForm] = solve_model_for_online_filter(true, xparam1, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

order = options_.order;
if order == 3 && ~ReducedForm.use_k_order_solver && ~isfield(ReducedForm, 'ghs3')
    ReducedForm.ghs3 = zeros(size(ReducedForm.steadystate));
end
mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
number_of_particles = online_opt.particles;
number_of_parameters = size(xparam1,1);
Y = dataset_.data;
sample_size = size(Y,1);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);
liu_west_delta = online_opt.liu_west_delta;

% Get initial conditions for the state particles
StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),ReducedForm.StateVectorMean);
if pruning && ~(order==1)
    % span initial augmented state vectors
    if order == 2
        StateVectors_ = StateVectors;
    elseif order == 3
        StateVectors_ = repmat(StateVectors,3,1);
    else
        error('Pruning is not available for orders > 3');
    end
else
    StateVectors_=[];
end

% parameters for the Liu & West filter
small_a = (3*liu_west_delta-1)/(2*liu_west_delta);
b_square = 1-small_a*small_a;

% Initialization of parameter particles
xparam = zeros(number_of_parameters,number_of_particles);
Prior = dprior(bayestopt_, options_.prior_trunc);
for i=1:number_of_particles
    info = 12042009;
    while info(1)
        candidate = Prior.draw();
        [info] = solve_model_for_online_filter(false, candidate, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        if ~info(1)
            xparam(:,i) = candidate(:);
        end
    end
end

% Initialization of the weights of particles.
weights = ones(1,number_of_particles)/number_of_particles;

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables;
mean_xparam = zeros(number_of_parameters,sample_size);
median_xparam = zeros(number_of_parameters,sample_size);
std_xparam = zeros(number_of_parameters,sample_size);
lb95_xparam = zeros(number_of_parameters,sample_size);
ub95_xparam = zeros(number_of_parameters,sample_size);

%% The Online filter
for t=1:sample_size
    if t>1
        fprintf('\nSubsample with %s first observations.\n\n', int2str(t))
    else
        fprintf('\nSubsample with only the first observation.\n\n')
    end
    % Moments of parameters particles distribution
    m_bar = xparam*(weights');
    temp = bsxfun(@minus,xparam,m_bar);
    sigma_bar = (bsxfun(@times,weights,temp))*(temp');
    if variance_update
        chol_sigma_bar = chol(b_square*sigma_bar)';
    end
    % Prediction (without shocks)
    fore_xparam = bsxfun(@plus,(1-small_a).*m_bar,small_a.*xparam);
    tau_tilde = zeros(1,number_of_particles);
    for i=1:number_of_particles
        % model resolution
        [info, M_, ReducedForm] = ...
            solve_model_for_online_filter(false, fore_xparam(:,i), options_, M_, estim_params_, bayestopt_, bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        if ~info(1)
            if ~ReducedForm.use_k_order_solver
                if pruning
                    if order == 2
                        state_variables_steady_state_ = ReducedForm.state_variables_steady_state;
                    elseif order == 3
                        state_variables_steady_state_ = repmat(ReducedForm.state_variables_steady_state,3,1);
                    else
                        error('Pruning is not available for orders > 3');
                    end
                end
            end
            % particle likelihood contribution
            yhat = bsxfun(@minus, StateVectors(:,i), ReducedForm.state_variables_steady_state);
            if ReducedForm.use_k_order_solver
                tmp = local_state_space_iteration_k(yhat, zeros(number_of_structural_innovations, 1), ReducedForm.dr, M_, options_, ReducedForm.udr);
            else
                if pruning
                    yhat_ = bsxfun(@minus,StateVectors_(:,i),state_variables_steady_state_);
                    if order == 2
                        [tmp,~] = local_state_space_iteration_2(yhat, zeros(number_of_structural_innovations, 1), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, yhat_, ReducedForm.steadystate, options_.threads.local_state_space_iteration_2);
                    elseif order == 3
                        [tmp,~] = local_state_space_iteration_3(yhat_, zeros(number_of_structural_innovations, 1), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
                    else
                    error('Pruning is not available for orders > 3');
                    end
                else
                    if order == 1
                        tmp = bsxfun(@plus,ReducedForm.constant,ReducedForm.ghx*yhat);
                    elseif order == 2
                        tmp = local_state_space_iteration_2(yhat, zeros(number_of_structural_innovations, 1), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, options_.threads.local_state_space_iteration_2);
                    elseif order == 3
                        tmp = local_state_space_iteration_3(yhat, zeros(number_of_structural_innovations, 1), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
                    else
                        error('Order > 3: use_k_order_solver should be set to true');
                    end
                end
            end
            % end
            PredictionError = bsxfun(@minus,Y(t,:)', tmp(mf1,:));
            % Replace Gaussian density with a Student density with 3 degrees of freedom for fat tails.
            z = sum(PredictionError.*(ReducedForm.H\PredictionError), 1) ;
            ddl = 3 ;
            tau_tilde(i) = weights(i)*(exp(gammaln((ddl + 1) / 2) - gammaln(ddl/2))/(sqrt(ddl*pi)*(1 + (z^2)/ddl)^((ddl + 1)/2))+1e-99) ;
        else
            tau_tilde(i) = 0 ;
        end
    end
    % particles selection
    tau_tilde = tau_tilde/sum(tau_tilde);
    indx = kitagawa(tau_tilde');
    StateVectors = StateVectors(:,indx);
    xparam = fore_xparam(:,indx);
    if pruning
        StateVectors_ = StateVectors_(:,indx);
    end
    w_stage1 = weights(indx)./tau_tilde(indx);
    % draw in the new distributions
    wtilde = zeros(1, number_of_particles);
    for i=1:number_of_particles
        info = 12042009;
        counter=0;
        while info(1) && counter <online_opt.liu_west_max_resampling_tries
            counter=counter+1;
            candidate = xparam(:,i) + chol_sigma_bar*randn(number_of_parameters, 1);
            if all(candidate>=bounds.lb) && all(candidate<=bounds.ub)
                % model resolution for new parameters particles
                [info, M_, ReducedForm] = ...
                    solve_model_for_online_filter(false, candidate, options_, M_, estim_params_, bayestopt_, bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
                if ~info(1)
                    xparam(:,i) = candidate ;
                    if order == 3 && ~ReducedForm.use_k_order_solver && ~isfield(ReducedForm, 'ghs3')
                        ReducedForm.ghs3 = zeros(size(ReducedForm.steadystate));
                    end
                    if ~ReducedForm.use_k_order_solver
                        if pruning
                            if order == 2
                                state_variables_steady_state_ = ReducedForm.state_variables_steady_state;
                                mf0_ = mf0;
                            elseif order == 3
                                state_variables_steady_state_ = repmat(ReducedForm.state_variables_steady_state,3,1);
                                mf0_ = repmat(mf0,1,3);
                                mask2 = number_of_state_variables+1:2*number_of_state_variables;
                                mask3 = 2*number_of_state_variables+1:number_of_state_variables;
                                mf0_(mask2) = mf0_(mask2)+size(ReducedForm.ghx,1);
                                mf0_(mask3) = mf0_(mask3)+2*size(ReducedForm.ghx,1);
                            else
                                error('Pruning is not available for orders > 3');
                            end
                        end
                    end
                    % Get covariance matrices and structural shocks
                    epsilon = chol(ReducedForm.Q)'*randn(number_of_structural_innovations, 1);
                    % compute particles likelihood contribution
                    yhat = bsxfun(@minus,StateVectors(:,i), ReducedForm.state_variables_steady_state);
                    if ReducedForm.use_k_order_solver
                        tmp = local_state_space_iteration_k(yhat, epsilon, ReducedForm.dr, M_, options_, ReducedForm.udr);
                    else
                        if pruning
                            yhat_ = bsxfun(@minus,StateVectors_(:,i), state_variables_steady_state_);
                            if order <= 2
                                [tmp, tmp_] = local_state_space_iteration_2(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, yhat_, ReducedForm.steadystate, options_.threads.local_state_space_iteration_2);
                            elseif order == 3
                                [tmp, tmp_] = local_state_space_iteration_3(yhat_, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ...
                                    ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
                            else
                                error('Pruning is not available for orders > 3');
                            end
                        else
                            if order == 1
                                tmp = bsxfun(@plus,ReducedForm.constant,ReducedForm.ghx*yhat)+ReducedForm.ghu*epsilon;
                            elseif order == 2
                                tmp = local_state_space_iteration_2(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, options_.threads.local_state_space_iteration_2);
                            elseif order == 3
                                tmp = local_state_space_iteration_3(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
                            else
                                error('Order > 3: use_k_order_solver should be set to true');
                            end
                        end
                    end
                    StateVectors(:,i) = tmp(mf0,:);
                    PredictionError = bsxfun(@minus,Y(t,:)', tmp(mf1,:));
                    wtilde(i) = w_stage1(i)*exp(-.5*(const_lik+log(det(ReducedForm.H))+sum(PredictionError.*(ReducedForm.H\PredictionError), 1)));
                end
            end
            if counter==online_opt.liu_west_max_resampling_tries
                fprintf('\nLiu & West online filter: I haven''t been able to solve the model in %u tries.\n',online_opt.liu_west_max_resampling_tries)
                fprintf('Liu & West online filter: The last error message was: %s\n',get_error_message(info))
                fprintf('Liu & West online filter: You can try to increase liu_west_max_resampling_tries, but most\n')
                fprintf('Liu & West online filter: likely there is an issue with the model.\n')
                error('Liu & West online filter: unable to solve the model.')
            end
        end
    end
    % normalization
    weights = wtilde/sum(wtilde);
    if variance_update && (neff(weights)<options_.particle.resampling.threshold*sample_size)
        variance_update = false;
    end
    % final resampling (not advised)
    if online_opt.second_resampling
        indx = kitagawa(weights);
        StateVectors = StateVectors(:,indx) ;
        if pruning
            StateVectors_ = StateVectors_(:,indx);
        end
        xparam = xparam(:,indx);
        weights = ones(1, number_of_particles)/number_of_particles;
        mean_xparam(:,t) = mean(xparam, 2);
        mat_var_cov = bsxfun(@minus, xparam, mean_xparam(:,t));
        mat_var_cov = (mat_var_cov*mat_var_cov')/(number_of_particles-1);
        std_xparam(:,t) = sqrt(diag(mat_var_cov));
        for i=1:number_of_parameters
            temp = sortrows(xparam(i,:)');
            lb95_xparam(i,t) = temp(0.025*number_of_particles);
            median_xparam(i,t) = temp(0.5*number_of_particles);
            ub95_xparam(i,t) = temp(0.975*number_of_particles);
        end
        if t==sample_size
            param = xparam ;
            save(sprintf('%s%sparameters_particles_final.mat', SimulationFolder, filesep()), 'param');
        end
    else
        mean_xparam(:,t) = xparam*(weights');
        mat_var_cov = bsxfun(@minus, xparam,mean_xparam(:,t));
        mat_var_cov = mat_var_cov*(bsxfun(@times, mat_var_cov, weights)');
        std_xparam(:,t) = sqrt(diag(mat_var_cov));
        for i=1:number_of_parameters
            temp = sortrows([xparam(i,:)' weights'], 1);
            cumulated_weights = cumsum(temp(:,2));
            pass1 = false;
            pass2 = false;
            pass3 = false;
            for j=1:number_of_particles
                if ~pass1 && cumulated_weights(j)>=0.025
                    lb95_xparam(i,t) = temp(j,1);
                    pass1 = true;
                end
                if  ~pass2 && cumulated_weights(j)>=0.5
                    median_xparam(i,t) = temp(j,1);
                    pass2 = true;
                end
                if ~pass3 && cumulated_weights(j)>=0.975
                    ub95_xparam(i,t) = temp(j,1);
                    pass3 = true;
                end
            end
        end
        if t==sample_size
            param = xparam(:,kitagawa(weights));
            save(sprintf('%s%sparameters_particles_final.mat', SimulationFolder, filesep()), 'param');
        end
    end
    save(sprintf('%s%sparameters_particles-%u.mat', SimulationFolder, filesep(), t), 'xparam');
    dyntable(options_,'', {'Parameter'; 'Lower Bound (95%)';'Mean';'Upper Bound (95%)'}, bayestopt_.name, [lb95_xparam(:,t), mean_xparam(:,t), ub95_xparam(:,t)], cellofchararraymaxlength(bayestopt_.name)+2, 10, 6)
    disp('')

end

%% Plot parameters trajectory
TeX = options_.TeX;

nr = ceil(sqrt(number_of_parameters)) ;
nc = floor(sqrt(number_of_parameters));
nbplt = 1 ;

if TeX
    fidTeX = fopen([M_.dname filesep 'online' filesep M_.fname '_param_traj.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by online_auxiliary_filter.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

for plt = 1:nbplt
    hh_fig = dyn_figure(options_.nodisplay,'Name','Parameters Trajectories');
    for k=1:number_of_parameters
        subplot(nr,nc,k)
        [name,texname] = get_the_name(k,TeX,M_,estim_params_,options_.varobs);
        % Draw the surface for an interval containing 95% of the particles.
        area(1:sample_size, ub95_xparam(k,:), 'FaceColor', [.9 .9 .9], 'BaseValue', min(lb95_xparam(k,:)));
        hold on
        area(1:sample_size, lb95_xparam(k,:), 'FaceColor', [1 1 1], 'BaseValue', min(lb95_xparam(k,:)));
        % Draw the mean of particles.
        plot(1:sample_size, mean_xparam(k,:), '-k', 'linewidth', 2)
        if TeX
            title(texname,'interpreter','latex')
        else
            title(name,'interpreter','none')
        end
        hold off
        axis tight
        drawnow
    end
    dyn_saveas(hh_fig, [M_.dname filesep 'online' filesep M_.fname '_param_traj' int2str(plt)], options_.nodisplay, options_.graph_format);
    if TeX
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s/online/%s_ParamTraj%s}\n',M_.dname,M_.fname,int2str(plt));
        fprintf(fidTeX,'\\caption{Parameters trajectories.}');
        fprintf(fidTeX,'\\label{Fig:ParametersPlots:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end
