function test_solve_with_custom_impulse_responses(M_, options_, oo_, x_hat_dash, G, J, F, curlyDs, curlyYs, varargin)
% Test heterogeneity.solve with custom impulse responses and SSJ validation
%
% INPUTS
% - M_ [structure]: Dynare model structure
% - options_ [structure]: Dynare options structure
% - oo_ [structure]: Dynare results structure
% - x_hat_dash [array]: custom impulse responses
% - G [structure]: SSJ general equilibrium Jacobians
% - J [structure]: SSJ heterogeneous-agent Jacobians
% - F [structure]: SSJ fake news matrices
% - curlyDs [structure]: SSJ aggregate impulse responses
% - curlyYs [structure]: SSJ heterogeneous impulse responses
% - varargin: Model-specific interpolation parameters and optional skip_init flag
%   For one-asset models: a_i, a_pi
%   For two-asset models: a_i, b_i, a_pi, b_pi

    % Set up interpolation weights based on model type
    if nargin == 11
        % One-asset model
        n_interp_args = 2;
        interpolation_indices = varargin{1};
        interpolation_weights = varargin{2};
        oo_.heterogeneity.mat.d.ind.a = interpolation_indices+1;
        oo_.heterogeneity.mat.d.w.a = interpolation_weights;
    elseif nargin == 13
        n_interp_args = 4;
        % Two-asset model
        a_i = varargin{1};
        b_i = varargin{2};
        a_pi = varargin{3};
        b_pi = varargin{4};
        oo_.heterogeneity.mat.d.ind.a = a_i+1;
        oo_.heterogeneity.mat.d.ind.b = b_i+1;
        oo_.heterogeneity.mat.d.w.a = a_pi;
        oo_.heterogeneity.mat.d.w.b = b_pi;
    end

    % Test solve function with custom impulse responses
    options_solve = options_.heterogeneity.solve;
    T = options_solve.truncation_horizon;
    oo_.heterogeneity.dr = heterogeneity.solve(M_, options_solve, oo_.heterogeneity, x_hat_dash);

    % Check that oo_.heterogeneity.dr structure contains expected fields
    assert(isfield(oo_.heterogeneity, 'dr'), 'oo_.heterogeneity.dr field missing');
    dr = oo_.heterogeneity.dr;
    assert(isfield(dr, 'G'), 'dr.G field missing');
    assert(isfield(dr, 'J_ha'), 'dr.J_ha field missing');
    assert(isfield(dr, 'F'), 'dr.F field missing');
    assert(isfield(dr, 'J'), 'dr.J field missing');
    assert(isfield(dr, 'curlyYs'), 'dr.curlyYs field missing');
    assert(isfield(dr, 'curlyDs'), 'dr.curlyDs field missing');

    % Test that results with custom impulse responses match SSJ reference
    tol = 1e-12;

    % curlyDs comparison - loop over all aggregate variables
    agg = fieldnames(curlyDs);
    for i=1:numel(agg)
        var = agg{i};
        M = abs(curlyDs.(var) - dr.curlyDs.(var));
        fprintf('curlyDs.%s vs SSJ reference residual: %.2e\n', var, max(M(:)));
        assert(max(M(:)) <= tol, sprintf('curlyDs.%s differs from SSJ reference with custom impulse responses', var));
    end

    % curlyYs comparison - loop over all variables and heterogeneous states
    het = fieldnames(dr.curlyYs.(agg{1}));
    for i=1:numel(agg)
        var = agg{i};
        for j=1:numel(het)
            x = het{j};
            M = abs(curlyYs.(var).(x).' - dr.curlyYs.(var).(x));
            fprintf('curlyYs.%s.%s vs SSJ reference residual: %.2e\n', var, x, max(M(:)));
            assert(max(M(:)) <= tol, sprintf('curlyYs.%s.%s differs from SSJ reference with custom impulse responses', var, x));
        end
    end

    % F matrices comparison - loop over heterogeneous states and aggregate variables
    J_fields = fieldnames(J);
    for j=1:numel(het)
        x = het{j};
        sum_x = sprintf('SUM_%s', x);
        for i=1:numel(agg)
            var = agg{i};
            % Map from SSJ naming (capital letters for J) to our naming
            if j <= numel(J_fields)
                M = abs(F.(x).(var) - dr.F.(sum_x).(var));
                fprintf('F.%s.%s vs SSJ reference residual: %.2e\n', sum_x, var, max(M(:)));
                assert(max(M(:)) <= T*tol, sprintf('F.%s.%s differs from SSJ reference with custom impulse responses', sum_x, var));
            end
        end
    end

    % Heterogeneous-agent Jacobians comparison
    for j=1:numel(het)
        x = het{j};
        sum_x = sprintf('SUM_%s', x);
        for i=1:numel(agg)
            var = agg{i};
            if j <= numel(J_fields)
                M = abs(J.(x).(var) - dr.J_ha.(sum_x).(var));
                fprintf('J_ha.%s.%s vs SSJ reference residual: %.2e\n', sum_x, var, max(M(:)));
                assert(max(M(:)) <= T*tol, sprintf('J_ha.%s.%s differs from SSJ reference with custom impulse responses', sum_x, var));
            end
        end
    end

    % General-equilibrium Jacobians comparison
    vars = intersect(fieldnames(G), fieldnames(dr.G));
    for i=1:numel(vars)
        v = vars{i};
        for j=1:M_.exo_nbr
            e = M_.exo_names{j};
            M = abs(dr.G.(v).(e) - G.(v).(e));
            fprintf('G.%s.%s vs SSJ reference residual: %.2e\n', v, e, max(M(:)));
            assert(max(M(:)) <= 1e-6, sprintf('G.%s.%s differs from SSJ reference', v, e));
        end
    end
end
