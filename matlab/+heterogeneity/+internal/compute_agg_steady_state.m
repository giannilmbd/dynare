function [y, yagg, x] = compute_agg_steady_state(M_, agg, x_bar_dash, Phi, hist, in_x)
   % Aggregate inputs:
   % - Aggregated heterogeneous variables
   Ix = x_bar_dash*Phi*hist;
   % - Extended aggregate endogenous variable vector at the steady state - %
   y = NaN(M_.endo_nbr,1);
   for i=1:M_.orig_endo_nbr
      var = M_.endo_names{i};
      y(i) = agg.(var); 
   end
   % - Variables resulting from an aggregation operator - %
   yagg = Ix(in_x);
   % - Exogenous variables - %
   x = zeros(M_.exo_nbr,1);
   % - Auxiliary variables - %
   if M_.set_auxiliary_variables
      fun_set_auxiliary_variables = str2func([M_.fname '.set_auxiliary_variables']);
      y = fun_set_auxiliary_variables(y, x, yagg, M_.params);
   end
end
