// --+ options: json=compute +--

// Declare the endogenous variables.
var ffr, unrate, cpi;

// Declare the exogenous variables.
varexo e_ffr, e_unrate, e_cpi;

parameters p_ffr_ffr_lag_1 p_ffr_ffr_lag_2 p_ffr_ffr_lag_3 p_ffr_unrate_lag_1 p_ffr_unrate_lag_2 p_ffr_unrate_lag_4 p_ffr_unrate_lag_5 p_ffr_cpi_lag_4 p_unrate_cpi_lag_1 p_unrate_cpi_lag_2 p_unrate_cpi_lag_3 p_unrate_cpi_lag_4 p_unrate_cpi_lag_5 p_unrate_cpi_lag_6 p_cpi_ffr_lag_1 p_cpi_ffr_lag_2 p_cpi_cpi_lag_2;

// Declare the model.
model(linear);

[name='ffr']
     ffr = p_ffr_ffr_lag_1*ffr_lag_1 + p_ffr_ffr_lag_2*ffr_lag_2 + p_ffr_ffr_lag_3*ffr_lag_3
           + p_ffr_unrate_lag_1*unrate_lag_1 + p_ffr_cpi_lag_4*cpi_lag_4 + e_ffr;

[name='unrate']
     unrate = unrate(-4)*p_ffr_unrate_lag_4 + unrate(-2)*p_ffr_unrate_lag_2 + unrate(-5)*p_ffr_unrate_lag_5 + cpi(-1)*p_unrate_cpi_lag_1 + cpi(-2)*p_unrate_cpi_lag_2 + cpi(-3)*p_unrate_cpi_lag_3 + cpi(-4)*p_unrate_cpi_lag_4 + cpi(-5)*p_unrate_cpi_lag_5 + cpi(-6)*p_unrate_cpi_lag_6 + e_unrate;

[name='cpi']
     cpi = ffr(-1)*p_cpi_ffr_lag_1 + ffr(-2)*p_cpi_ffr_lag_2 + cpi(-2)*p_cpi_cpi_lag_2 + e_cpi;

end;

// Set true values for the parameters (randomly, I write directly in M_.params
// because I am too lazy to write all the parameter names).
M_.params = .1*randn(M_.param_nbr, 1);

shocks;
var e_ffr; stderr .01;
var e_unrate; stderr .01;
var e_cpi; stderr .01;
end;

// Create a dataset by simulating the model.
oo_ = simul_backward_model([], 1000, options_, M_, oo_);

// Put all the simulated data in a dseries object.
ds1 = dseries(transpose(oo_.endo_simul), 1900Q1, M_.endo_names);


// Select a subsample for estimation
ds1 = ds1(1951Q1:1990Q1);


// Print true values of the parameters (for comparison with the estimation results).
dyn_table('True values for the parameters (DGP)', {}, {}, cellstr(M_.param_names), ...
    {'Parameter Value'}, 4, M_.params);


/* Run the estimation of equation ffr and display the results. Note that the corresponding
** values in M_.params will be updated.
**
** estimate command
** ----------------
**
** First argument is used to specify the estimation method
** Second argument is used to specify the dataset (an existing dseries object).
** Third argument (and following) is(are) the name(s) of equation(s) to be estimated.
*/

estimate method(ols) data(ds1) ffr;
