profile off;
if ~isoctave
profData = profile('info');
% profile viewer % this opens the GUI

% define the functions to look for that call 'dsge_likelihood.m' during optimization
switch options_.mode_compute
    case {1, 3, 7, 12, 102}
        % These are optimization algorithms provided by MATLAB
        % The direct caller of our function is not documented and can change
        % across MATLAB versions.
        % In any case, we do not have to test that MATLAB optimizers correctly
        % counts function calls, that’s not the purpose of our testsuite.
        fun = {};
    case 2
        fun = {'simulated_annealing'};
    case 4
        fun = {'penalty_objective_function'};
    case 5
        fun = {'penalty_objective_function'};
    case 6
        fun = {'gmhmaxlik','gmhmaxlik_core'};
    case 8
        fun = {'simplex_optimization_routine', 'simplex_optimization_routine>simplex_initialization'};
    case 9
        fun = {'cmaes'};
    case 10
        fun = {'simpsa>CALCULATE_COST'};
    case 101
        if options_.analytic_derivation
            fun = {'analytic_gradient_wrapper'};
        else
            fun = {'solvopt', 'apprgrdn'};
        end
    otherwise
        fun = {'penalty_objective_function'};
end
% find dsge_likelihood function
dsge_idx = find(strcmp({profData.FunctionTable.FunctionName}, 'dsge_likelihood'));
dsge_func = profData.FunctionTable(dsge_idx);
% get all parent function names using their indices
parent_indices = [dsge_func.Parents.Index];
parent_names = {profData.FunctionTable(parent_indices).FunctionName};
% find matching functions and sum their NumCalls
matching_idx = ismember(parent_names, fun);
NumCalls = sum([dsge_func.Parents(matching_idx).NumCalls]);

if options_.analytic_derivation
    fprintf('OPTIMIZATION INFO WITH ANALYTICAL DERIVATIVES\n');
else
    fprintf('OPTIMIZATION INFO\n');
end
fprintf('- funcCount: %d\n', oo_.posterior.optimization.optimization_info.funcCount);
fprintf('- runtime: %s\n', dynsec2hms(oo_.posterior.optimization.optimization_info.runtime));
fprintf('- iterations: %d\n', oo_.posterior.optimization.optimization_info.iterations);
fprintf('- exitflag: %d\n', oo_.posterior.optimization.optimization_info.exitflag);
fprintf('- termination message: %s\n', oo_.posterior.optimization.optimization_info.message);

if ~isempty(fun) && ~isequal(oo_.posterior.optimization.optimization_info.funcCount,NumCalls)
    error('counting of function evaluations is not correct');
end
end
