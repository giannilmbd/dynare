% Tests for preprocessor capabilities of skew in estimated_params blocks
% Note: manual tests for error messages by uncommenting specific lines
% -------------------------------------------------------------------------

% Copyright © 2025 Dynare Team
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

@#include "ireland2004_model.inc"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: check estim_params with concatenating %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimated_params;
skew eta_a, normal_pdf, 1, 0.10, 0.11, 1.11, 111;
skew eta_e, 2.22, 2.00, 22.22, beta_pdf, 2, 0.22;
skew eta_r, gamma_pdf, 3, 0.30;
%skew yhat, inv_gamma_pdf, 5, 0.50; % triggers error because variable cannot be endogenous (yhat)
%skew eta_u, normal_pdf, 38, 0.38; % triggers error because of undeclared eta_u
%skew ALPHA_PI, normal_pdf, 38, 0.38; % triggers error because variable cannot be parameter (ALPHA_PI)
end;

fprintf('\nPreprocessing of estimated_params block...\n')
fprintf('%s estim_params_.skew_exo %s\n\n', repmat('=',1,40), repmat('=',1,40));
disp(estim_params_.skew_exo)
fprintf('%s\n', repmat('=',1,104));

skew_exo = [
    1     nan    -Inf    Inf    3    1    0.10   0.11  1.11  111;
    2    2.22       2  22.22    1    2    0.22   nan   nan   nan;
    4     nan    -Inf    Inf    2    3    0.30   nan   nan   nan;
];
if ~isequaln(skew_exo, estim_params_.skew_exo)
    error('estim_params block is preprocessed incorrectly');
else
    fprintf('...successful!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: check estimated_params_init %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params_init;
skew eta_a, 555;
skew eta_e, 222;
skew eta_z, 666; % displays message: The skewness of eta_z is not estimated (the value provided in estimated_params_init is not used).
end;

skew_exo(1, 2) = 555;
skew_exo(2, 2) = 222;

fprintf('\nPreprocessing of estimated_params_init block...\n')
fprintf('%s estim_params_.skew_exo %s\n\n', repmat('=',1,40), repmat('=',1,40));
disp(estim_params_.skew_exo)
fprintf('%s\n', repmat('=',1,104));
if ~isequaln(skew_exo, estim_params_.skew_exo)
    error('estimated_params_init block is preprocessed incorrectly');
else
    fprintf('...successful!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: check estimated_params_bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params_bounds;
skew eta_a, 5, 55;
skew eta_e, 2, 22;
end;

skew_exo(1, [3 4]) = [5 55];
skew_exo(2, [3 4]) = [2 22];

fprintf('\nPreprocessing of estimated_params_bounds block...\n')
fprintf('%s estim_params_.skew_exo %s\n\n', repmat('=',1,40), repmat('=',1,40));
disp(estim_params_.skew_exo)
fprintf('%s\n', repmat('=',1,104));
if ~isequaln(skew_exo, estim_params_.skew_exo)
    error('estimated_params_bounds block is preprocessed incorrectly');
else
    fprintf('...successful!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: check estimated_params_remove %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params_remove;
skew eta_e;
end;

skew_exo(2,:) = [];

fprintf('\nPreprocessing of estimated_params_remove block...\n')
fprintf('%s estim_params_.skew_exo %s\n\n', repmat('=',1,40), repmat('=',1,40));
disp(estim_params_.skew_exo)
fprintf('%s\n', repmat('=',1,104));
if ~isequaln(skew_exo, estim_params_.skew_exo)
    error('estimated_params_remove block is preprocessed incorrectly');
else
    fprintf('...successful!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 5: check estimated_params(overwrite) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params(overwrite);
skew eta_z, inv_gamma_pdf, 4, 0.40;
end;
skew_exo = [3, nan, -Inf, Inf, 4, 4, 0.40, nan, nan, nan];

fprintf('\nPreprocessing of estim_params(overwrite) block...\n')
fprintf('%s estim_params_.skew_exo %s\n\n', repmat('=',1,40), repmat('=',1,40));
disp(estim_params_.skew_exo)
fprintf('%s\n', repmat('=',1,104));
if ~isequaln(skew_exo, estim_params_.skew_exo)
    error('estim_params(overwrite) block is preprocessed incorrectly');
else
    fprintf('...successful!\n')
end