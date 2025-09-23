function test_solve_functionality(M_, oo_, G, varargin)
% Test heterogeneity.solve functionality and SSJ validation
%
% INPUTS
% - M_ [structure]: Dynare model structure
% - options_ [structure]: Dynare options structure
% - oo_ [structure]: Dynare results structure
% - ss [structure]: steady-state structure
% - G [structure]: SSJ general equilibrium Jacobians for validation
% - varargin [optional]:
%     varargin{1} - flag: if true, assert on validation errors (default: true)
    if nargin>3 && ~isempty(varargin{1})
        flag = varargin{1};
    else
        flag = true;
    end

    % Check that oo_.heterogeneity.dr structure contains expected fields
    assert(isfield(oo_.heterogeneity, 'dr'), 'oo_.heterogeneity.dr field missing');
    dr = oo_.heterogeneity.dr;
    assert(isfield(dr, 'G'), 'dr.G field missing');
    assert(isfield(dr, 'J_ha'), 'dr.J_ha field missing');
    assert(isfield(dr, 'F'), 'dr.F field missing');
    assert(isfield(dr, 'J'), 'dr.J field missing');
    assert(isfield(dr, 'curlyYs'), 'dr.curlyYs field missing');
    assert(isfield(dr, 'curlyDs'), 'dr.curlyDs field missing');

    % Test comparison with SSJ results
    % General-equilibrium Jacobians comparison
    vars = intersect(fieldnames(G), fieldnames(dr.G));
    for i=1:numel(vars)
        v = vars{i};
        % Loop over shocks
        for j=1:M_.exo_nbr
            e = M_.exo_names{j};
            M = abs(dr.G.(v).(e) - G.(v).(e));
            fprintf('G.%s.%s vs SSJ reference residual: %.2e\n', v, e, max(M(:)));
            if flag
                assert(max(M(:)) <= 1e-3, sprintf('G.%s.%s differs from SSJ reference', v, e));
            end
        end
    end
end
