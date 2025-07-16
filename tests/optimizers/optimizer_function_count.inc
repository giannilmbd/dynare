profile off; profData = profile('info');
% profile viewer % this opens the GUI

if options_.mode_compute == 12
    % number of calls of penalty_objective_function.m from profiler
    NumCalls = profData.FunctionTable(find(strcmp({profData.FunctionTable.FunctionName}, 'penalty_objective_function'))).NumCalls;
else
    % number of calls of dsge_likelihood.m from profiler
    NumCalls = profData.FunctionTable(find(strcmp({profData.FunctionTable.FunctionName}, 'dsge_likelihood'))).NumCalls;
    % common calling functions:
    % -  (1 or 2 with analytic_derivation)
    % - 
    verbatim;
    offset = 1; % evaluation in initial_estimation_checks.m
    if (isnumeric(options_.mode_compute) && ~ismember(options_.mode_compute, [9, 10, 102])) || (ischar(options_.mode_compute) && strcmp(options_.mode_compute, 'optimizer_function_wrapper'))
        offset = offset + 1; % evaluation in dynare_estimation_1.m
    end
    if options_.analytic_derivation
        offset = offset + 1; % additional evaluation in dynare_estimation_1.m
    else
        if options_.mode_compute ~= 6
            n = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn + estim_params_.np;
            offset = offset + 1 + 2*n + n*(n-1); % evaluations in hessian.m
        end
    end
    end;
    NumCalls = NumCalls - offset;
    
    % number of calls of dsge_likelihood.m from optimizers
    if ~isequal(oo_.posterior.optimization.optimization_info.funcCount,NumCalls)
        error('counting of function evaluations of dsge_likelihood.m is not correct');
    else
        if options_.analytic_derivation
            fprintf('OPTIMIZATION INFO WITH ANALYTICAL DERIVATIVES\n');
        else
            fprintf('OPTIMIZATION INFO\n');
        end
        fprintf('dsge_likelihood was called %d times by the optimizer with a runtime of %s\n', NumCalls, dynsec2hms(oo_.posterior.optimization.optimization_info.runtime));
    end
end