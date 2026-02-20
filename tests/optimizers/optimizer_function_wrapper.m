function [opt_par_values,fval,exitflag,output]=optimizer_function_wrapper(objective_function_handle,start_par_value,varargin)
% function [opt_par_values,fval,exitflag,output]=optimizer_function_wrapper(objective_function_handle,start_par_value,varargin)
% Demonstrates how to invoke external optimizer for mode_computation

%set options of optimizer
H0 = 1e-4*eye(length(start_par_value),length(start_par_value));
nit=1000;
crit = 1e-7;
numgrad = 2;
epsilon = 1e-6;
analytic_grad=[];
Verbose=1;
Save_files=1;
%call optimizer
[fval,opt_par_values,grad,hessian_mat,output.iterations,output.funcCount,exitflag,output.message] = ...
    csminwel1(objective_function_handle, start_par_value, H0, analytic_grad, crit, nit, numgrad, epsilon, Verbose,Save_files, varargin{:});
end