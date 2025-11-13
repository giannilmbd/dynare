function a=set_Kalman_starting_values(a,M_,dr,options_,bayestopt_)
% function a=set_Kalman_starting_values(a,M_,dr,options_,bayestopt_)
% Sets initial states guess for Kalman filter/smoother based on M_.filter_initial_state
%
% INPUTS
%   o a             [double]   (p*1) vector of states
%   o M_            [structure] describing the model
%   o dr            [structure] storing the decision rules
%   o options_      [structure] describing the options
%   o bayestopt_    [structure] describing the priors
%
% OUTPUTS
%   o a             [double]    (p*1) vector of set initial states

if isfield(M_,'filter_initial_state') && ~isempty(M_.filter_initial_state)
    state_indices=dr.order_var(dr.restrict_var_list(bayestopt_.mf0));
    for ii=1:size(state_indices,1)
        if ~isempty(M_.filter_initial_state{state_indices(ii),1})
            if options_.loglinear && ~options_.logged_steady_state
                a(bayestopt_.mf0(ii)) = log(eval(M_.filter_initial_state{state_indices(ii),2})) - log(dr.ys(state_indices(ii)));
            elseif ~options_.loglinear && ~options_.logged_steady_state
                a(bayestopt_.mf0(ii)) = eval(M_.filter_initial_state{state_indices(ii),2}) - dr.ys(state_indices(ii));
            else
                error('The steady state is logged. This should not happen. Please contact the developers')
            end
        end
    end
end
