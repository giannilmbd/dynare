function [x_t_plus_1, x_t_plus_1_pruned_state_space]=iterate_law_of_motion(initial_state,epsilon,ReducedForm,M_,options_,use_k_order_solver,pruning,initial_pruned_state_space)
% function x_t_plus_1=iterate_law_of_motion(initial_state_vector,epsilon,ReducedForm,M_,options_,use_k_order_solver,pruning,initial_pruned_state_space)
% Inputs
%  - initial_state              [double]    initial states vector
%  - epsilon                    [double]    shock vector
%  - ReducedForm                [structure] decision rules
%  - M_                         [structure] describing the model
%  - options_                   [structure] describing the options
%  - use_k_order_solver         [boolean]   indicator whether the k_order_solver should be used
%  - pruning                    [boolean]   indicator whether pruning should be used
%  - initial_pruned_state_space [double]    components of the augmented pruned state space (required for pruning)
%
% Outputs:
% - x_t_plus_1                  [double]    simulated series next period
% - x_t_plus_1_pruned_state_space [double]  components of the pruned state space of the simulated series next period

if pruning
    if nargin<8
        error('iterate_law_of_motion: with pruning, the components of the augmented pruned state space are an required input')
    end
end

n_shock_vectors=size(epsilon,2);

x_t_plus_1_pruned_state_space=[]; %initialize to run for order=1 and pruning

yhat = bsxfun(@minus, initial_state, ReducedForm.state_variables_steady_state);
if size(initial_state,2)==1 && n_shock_vectors>1
    yhat=repmat(yhat,1,n_shock_vectors);
end
if use_k_order_solver
    x_t_plus_1 = local_state_space_iteration_k(yhat, zeros(number_of_structural_innovations, 1), ReducedForm.dr, M_, options_, ReducedForm.udr);
else
    if options_.order == 3 && ~isfield(ReducedForm, 'ghs3')
        ReducedForm.ghs3 = zeros(size(ReducedForm.steadystate));
    end
    if options_.order == 1
        x_t_plus_1 = bsxfun(@plus,ReducedForm.constant,ReducedForm.ghx*yhat)+ReducedForm.ghu*epsilon;
    else
        if pruning
            if options_.order == 2
                yhat_pruned_state_space = bsxfun(@minus,initial_pruned_state_space,ReducedForm.state_variables_steady_state);
                [x_t_plus_1, x_t_plus_1_pruned_state_space]= local_state_space_iteration_2(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ...
                    ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, yhat_pruned_state_space, ReducedForm.steadystate, ...
                    options_.threads.local_state_space_iteration_2);
            elseif options_.order == 3
                yhat_pruned_state_space = bsxfun(@minus,initial_pruned_state_space,repmat(ReducedForm.state_variables_steady_state,3,1));
                [x_t_plus_1, x_t_plus_1_pruned_state_space]= local_state_space_iteration_3(yhat_pruned_state_space, epsilon, ReducedForm.ghx, ReducedForm.ghu, ...
                    ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ...
                    ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ...
                    ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
            else
                error('Pruning is not available for orders > 3');
            end
        else
            if options_.order == 2
                x_t_plus_1 = local_state_space_iteration_2(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ...
                    ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, options_.threads.local_state_space_iteration_2);
            elseif options_.order == 3
                x_t_plus_1 = local_state_space_iteration_3(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ...
                    ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ...
                    ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ...
                    ReducedForm.steadystate, ReducedForm.ghs3, options_.threads.local_state_space_iteration_3, pruning);
            else
                error('Order > 3: use_k_order_solver should be set to true');
            end
        end
    end
end
