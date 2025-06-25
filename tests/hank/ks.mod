heterogeneity_dimension households;

var(heterogeneity=households) 
    c  // Consumption
    a  // Assets
    e  // Idiosyncratic (log-)efficiency
    Va // Derivative of the value function w.r.t assets
;

varexo(heterogeneity=households) eps_e; // Idiosyncratic efficiency shock

var
    r // Rate of return on capital net of depreciation
    w // Wage rate
    Y // Aggregate output
    K // Aggregate capital
    Z // Aggregate productivity
;
varexo eps_Z; // Aggregate productivity shock

parameters
    L     // Labor
    alpha // Share of capital in production fuction
    beta  // Subjective discount rate of houselholds
    delta // Capital depreciation rate
    eis   // Elasticity of intertemporal substitution
    rho_e // Earning shock persistence
    sig_e // Earning shock innovation std err
    rho_Z // Aggregate TFP shock persistence
    sig_Z // Aggregate TFP shock innovation std err
    Z_ss  // Aggregate TFP shock average value
;

model(heterogeneity=households);
    c^(-1/eis)-beta*Va(+1)=0 ⟂ a>=0;
    (1+r)*a(-1)+w*e-c-a; 
    Va = (1+r)*c^(-1/eis);
    log(e) - rho_e*log(e(-1)) - eps_e;
end;

model;
    Z * K(-1)^alpha * L^(1 - alpha) - Y;
    alpha * Z * (K(-1) / L)^(alpha - 1) - delta - r;
    (1 - alpha) * Z * (K(-1) / L)^alpha - w;
    K - SUM(a);
    log(Z) - rho_Z*log(Z(-1)) - (1-rho_Z)*log(Z_ss) - eps_Z;
end;

alpha = 0.11;
beta = 0.9819527880123727;
eis = 1;
delta = 0.025;
L = 1;
rho_e = 0.966;
sig_e = 0.5*sqrt(1-rho_e^2);
rho_Z = 0.8;
sig_Z = 0.014;
Z_ss = 0.8816;

shocks(heterogeneity=households);
    var eps_e; stderr sig_e;
end;

shocks;
    var eps_Z; stderr sig_Z;
end;

verbatim;
    skipline()
    disp('*** TESTING: hank.check_steady_state_input.m - Discretized AR(1) case ***');
    load('ks_ar1_ss.mat');
    verbose = true;
    options_.hank.nowarningredundant = false;
    options_.hank.nowarningdgrids = false;
    testFailed = 0;

    %% === OUTPUTS ===
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, base_struct);
        if ~(sizes.n_e == numel(fieldnames(base_struct.shocks.grids)) && ...
             sizes.n_a == numel(fieldnames(base_struct.pol.grids)) && ...
             sizes.N_e == numel(base_struct.shocks.grids.e) && ...
             sizes.n_pol == M_.heterogeneity(1).endo_nbr && ...
             sizes.agg == numel(fieldnames(base_struct.agg)) && ...
             sizes.shocks.e == numel(base_struct.shocks.grids.e) && ...
             sizes.pol.N_a == numel(base_struct.pol.grids.a) && ...
             sizes.pol.states.a == sizes.pol.N_a && ...
             sizes.d.N_a == sizes.pol.N_a && sizes.d.states.a == sizes.pol.states.a)
          testFailed = testFailed+1;
          dprintf('The `sizes` output is not correct!');
        end
        if ~(iscolumn(out_ss.shocks.grids.e) && iscolumn(out_ss.pol.grids.a) && ...
             iscolumn(out_ss.d.grids.a) && ...
             size(out_ss.shocks.grids.e,1) == sizes.shocks.e && ...
             size(out_ss.pol.grids.a,1) == sizes.pol.states.a && ...
             size(out_ss.d.grids.a,1) == sizes.d.states.a && ...
             numel(fieldnames(out_ss.agg)) == sizes.agg)
          testFailed = testFailed+1;
          dprintf('The `out_ss` output is not correct!');
        end
    catch ME
        testFailed = testFailed+2;
        dprintf('Outputs `sizes` or `out_ss` are not correct!')
    end

    %% === NON-STRUCT FIELDS ===

    % ss is not a struct
    ss = 123;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss is not a struct', testFailed, verbose);

    % ss.shocks is not a struct
    ss = base_struct;
    ss.shocks = 1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks is not a struct', testFailed, verbose);

    % ss.shocks.grids is not a struct
    ss = base_struct;
    ss.shocks.grids = 42;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks.grids is not a struct', testFailed, verbose);

    % ss.shocks.Pi is not a struct
    ss = base_struct;
    ss.shocks.Pi = "not_a_struct";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks.Pi is not a struct', testFailed, verbose);

    % ss.shocks.w is not a struct
    ss = base_struct;
    ss.shocks = rmfield(ss.shocks, 'Pi');
    ss.shocks.w = "not_a_struct";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks.w is not a struct', testFailed, verbose);

    % ss.pol is not a struct
    ss = base_struct;
    ss.pol = pi;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol is not a struct', testFailed, verbose);

    % ss.pol.grids is not a struct
    ss = base_struct;
    ss.pol.grids = "grid";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol.grids is not a struct', testFailed, verbose);

    % ss.pol.values is not a struct
    ss = base_struct;
    ss.pol.values = "values";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol.values is not a struct', testFailed, verbose);

    % ss.d is not a struct
    ss = base_struct;
    ss.d = 1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.d is not a struct', testFailed, verbose);

    % ss.d.grids is not a struct
    ss = base_struct;
    ss.d.grids = "grid";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.d.grids is not a struct', testFailed, verbose);

    % ss.agg is not a struct
    ss = base_struct;
    ss.agg = 0;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.agg is not a struct', testFailed, verbose);

    %% === SHOCKS TESTS ===

    % Missing `shocks` field
    ss = base_struct;
    ss = rmfield(ss, 'shocks');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `shocks` field', testFailed, verbose);

    % Missing `shocks.grids`
    ss = base_struct;
    ss.shocks = rmfield(ss.shocks, 'grids');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `shocks.grids`', testFailed, verbose);

    % Wrong type for `shocks.grids.e`
    ss = base_struct;
    ss.shocks.grids.e = ones(3,3);  % invalid: cell array
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'shocks.grids.e are not dense real vectors', testFailed, verbose);

    % Missing both `shocks.Pi` and `shocks.w`
    ss = base_struct;
    ss.shocks = rmfield(ss.shocks, 'Pi');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'No discretization method (Pi or w) provided', testFailed, verbose);

    % Both `shocks.Pi` and `shocks.w` provided (mutual exclusivity)
    ss = base_struct;
    ss.shocks.w = struct('e', ones(size(ss.shocks.Pi.e,1),1));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Both Pi and w should not coexist', testFailed, verbose);

    % AR(1) Tests
    % Shock grid not in var(heterogeneity=) list
    ss = base_struct;
    ss.shocks.grids.fake_shock = ss.shocks.grids.e;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Shock grid symbol not in heterogeneity declaration', testFailed, verbose);

    % Missing Markov matrix for a shock
    ss = base_struct;
    ss.shocks.Pi = rmfield(ss.shocks.Pi, 'e');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing Pi entry for a declared shock', testFailed, verbose);

    % Wrong type for `shocks.grids.Pi.e`
    ss = base_struct;
    ss.shocks.Pi.e = speye(numel(ss.shocks.grids.e));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'shocks.Pi.e are not dense real matrices', testFailed, verbose);

    % Markov matrix wrong row number
    ss = base_struct;
    ss.shocks.Pi.e = rand(3, 7);  % Should be square matching grid
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Pi row count mismatch with grid', testFailed, verbose);

    % Markov matrix wrong column number
    ss = base_struct;
    ss.shocks.Pi.e = rand(7, 3);  % Should be square matching grid
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Pi column count mismatch with grid', testFailed, verbose);

    % Markov matrix has negative elements
    ss = base_struct;
    ss.shocks.Pi.e(1,1) = -0.1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Negative elements in Pi', testFailed, verbose);

    % Markov matrix does not sum to 1
    ss = base_struct;
    ss.shocks.Pi.e = ss.shocks.Pi.e * 2;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Rows of Pi not summing to 1', testFailed, verbose);


    %% === POLICY TESTS ===

    % Missing `pol` field
    ss = base_struct;
    ss = rmfield(ss, 'pol');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol` field', testFailed, verbose);

    % Missing `pol.grids`
    ss = base_struct;
    ss.pol = rmfield(ss.pol, 'grids');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol.grids`', testFailed, verbose);

    % Missing one state grid
    ss = base_struct;
    ss.pol.grids = rmfield(ss.pol.grids, 'a');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing state grid for a', testFailed, verbose);

    % Wrong type for `pol.grids.a`
    ss = base_struct;
    ss.pol.grids.a = ones(3, 3);  % should be a vector
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.grids.a should be a vector', testFailed, verbose);

    % Missing `pol.values`
    ss = base_struct;
    ss.pol = rmfield(ss.pol, 'values');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol.values`', testFailed, verbose);

    % Missing policy value for declared variable
    ss = base_struct;
    ss.pol.values = rmfield(ss.pol.values, 'a');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing policy value for a', testFailed, verbose);

    % Wrong type for `pol.values.a`
    ss = base_struct;
    ss.pol.values.a = rand(2, 2, 2); 
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.values.a should be a matrix', testFailed, verbose);

    % Incompatible policy matrix row size (shocks × states)
    ss = base_struct;
    ss.pol.values.a = rand(99, size(ss.pol.values.a, 2));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of rows in policy matrix', testFailed, verbose);

    % Incompatible policy matrix column size
    ss = base_struct;
    ss.pol.values.a = rand(size(ss.pol.values.a, 1), 99);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of columns in policy matrix', testFailed, verbose);

    %% === DISTRIBUTION TESTS ===

    % Missing `d` field
    ss = base_struct;
    ss = rmfield(ss, 'd');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing distribution structure `d`', testFailed, verbose);

    % Wrong type for `ss.d.grids.a`
    ss = base_struct;
    ss.d.grids.a = rand(3, 3);  % not a vector
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.grids.a should be a vector', testFailed, verbose);

    % Missing `d.hist`
    ss = base_struct;
    ss.d = rmfield(ss.d, 'hist');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing d.hist', testFailed, verbose);

    % Wrong type for `ss.d.hist`
    ss = base_struct;
    ss.d.hist = rand(2,2,2);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.hist should be a matrix', testFailed, verbose);

    % Wrong row size in `d.hist`
    ss = base_struct;
    ss.d.hist = rand(99, size(ss.d.hist,2));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of rows in d.hist', testFailed, verbose);

    % Wrong column size in `d.hist`
    ss = base_struct;
    ss.d.hist = rand(size(ss.d.hist,1), 99);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of columns in d.hist', testFailed, verbose);

    %% === AGGREGATES TESTS ===

    % Missing `agg` field
    ss = base_struct;
    ss = rmfield(ss, 'agg');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `agg` field', testFailed, verbose);

    % Remove one required aggregate variable
    ss = base_struct;
    ss.agg = rmfield(ss.agg, 'r');  % remove interest rate
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing aggregate variable r', testFailed, verbose);

    % Wrong type for `ss.agg.r`
    ss = base_struct;
    ss.agg.r = 'non-numeric';
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'agg.r should be numeric', testFailed, verbose);

    %% === PERMUTATION TESTS ===

    % pol.shocks is not a string or cellstr
    ss = base_struct;
    ss.pol.shocks = 123;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.shocks is not a string array or cellstr', testFailed, verbose);

    % pol.shocks has wrong names
    ss = base_struct;
    ss.pol.shocks = {'wrong_name'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.shocks fieldnames is not consistent', testFailed, verbose);

    % pol.states is not a string or cellstr
    ss = base_struct;
    ss.pol.states = 42;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.states is not a string array or cellstr', testFailed, verbose);

    % pol.states has wrong names
    ss = base_struct;
    ss.pol.states = {'bad_state'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.states not matching state_var', testFailed, verbose);

    % d.shocks is not a string or cellstr
    ss = base_struct;
    ss.d.shocks = 3.14;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.shocks is not a string array or cellstr', testFailed, verbose);

    % d.shocks inconsistent with pol.shocks
    ss = base_struct;
    ss.d.shocks = {'ghost_shock'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.shocks not matching pol.shocks', testFailed, verbose);

    % d.states is not a string or cellstr
    ss = base_struct;
    ss.d.states = true;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.states is not a string array or cellstr', testFailed, verbose);

    % d.states inconsistent with pol.states
    ss = base_struct;
    ss.d.states = {'phantom_state'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.states not matching pol.states', testFailed, verbose);

    %% === REDUNDANCY TESTS ===

    % shocks.Pi contains redundant entry (not in shocks.grids)
    ss = base_struct;
    ss.shocks.Pi.extra = eye(length(ss.shocks.grids.e));  % not listed in shocks.grids
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant shocks.Pi entry accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant shocks.Pi entry: %s\n', ME.message);
    end

    % pol.values contains redundant field (not in model's pol_symbs)
    ss = base_struct;
    ss.pol.values.redundant_var = ones(size(ss.pol.values.a));
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant pol.values field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant pol.values field: %s\n', ME.message);
    end

    % d.grids contains redundant field (not a declared state)
    ss = base_struct;
    ss.d.grids.not_a_state = [0;1;2];  % not in M_.heterogeneity.state_var
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant d.grids field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant d.grids field: %s\n', ME.message);
    end

    % agg contains extra field (not in M_.endo_names)
    ss = base_struct;
    ss.agg.extra_agg = 999;
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant agg field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant agg field: %s\n', ME.message);
    end

    % Test the initialize_steady_state routine
    try
        oo_ = hank.initialize_steady_state(M_, options_, oo_, base_struct);
        disp('✔ Initialization of the steady state succeeded!');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from initialize_steady_state: %s\n', ME.message);
    end

    skipline()
    disp('*** TESTING: hank.check_steady_state_input.m - Discretized i.i.d case ***');
    load 'ks_iid_ss.mat';
    verbose = true;

    %% === OUTPUTS ===
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, base_struct);
        N_a = numel(base_struct.pol.grids.a)*numel(base_struct.pol.grids.e);
        if ~(sizes.n_e == numel(fieldnames(base_struct.shocks.grids)) && ...
             sizes.n_a == numel(fieldnames(base_struct.pol.grids)) && ...
             sizes.N_e == numel(base_struct.shocks.grids.eps_e) && ...
             sizes.n_pol == M_.heterogeneity(1).endo_nbr && ...
             sizes.agg == numel(fieldnames(base_struct.agg)) && ...
             sizes.shocks.eps_e == numel(base_struct.shocks.grids.eps_e) && ...
             sizes.pol.N_a == N_a && sizes.pol.states.a == numel(base_struct.pol.grids.a) && ...
             sizes.pol.states.e == numel(base_struct.pol.grids.e) && ...
             sizes.d.N_a == sizes.pol.N_a && sizes.d.states.a == sizes.pol.states.a)
          testFailed = testFailed+1;
          dprintf('The `sizes` output is not correct!');
        end
        if ~(iscolumn(out_ss.shocks.grids.eps_e) && iscolumn(out_ss.pol.grids.a) && ...
             iscolumn(out_ss.pol.grids.e) && iscolumn(out_ss.d.grids.a) && ...
             iscolumn(out_ss.d.grids.e) && ...
             size(out_ss.shocks.grids.eps_e,1) == sizes.shocks.eps_e && ...
             size(out_ss.pol.grids.a,1) == sizes.pol.states.a && ...
             size(out_ss.pol.grids.e,1) == sizes.pol.states.e && ...
             size(out_ss.d.grids.a,1) == sizes.d.states.a && ...
             size(out_ss.d.grids.e,1) == sizes.d.states.e && ...
             numel(fieldnames(out_ss.agg)) == sizes.agg)
          testFailed = testFailed+1;
          dprintf('The `out_ss` output is not correct!');
        end
    catch ME
        testFailed = testFailed+2;
        dprintf('Outputs `sizes` or `out_ss` are not correct!')
    end

    %% === NON-STRUCT FIELDS ===

    % ss is not a struct
    ss = 123;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss is not a struct', testFailed, verbose);

    % ss.shocks is not a struct
    ss = base_struct;
    ss.shocks = 1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks is not a struct', testFailed, verbose);

    % ss.shocks.grids is not a struct
    ss = base_struct;
    ss.shocks.grids = 42;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks.grids is not a struct', testFailed, verbose);

    % ss.shocks.w is not a struct
    ss = base_struct;
    ss.shocks.w = "not_a_struct";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.shocks.w is not a struct', testFailed, verbose);

    % ss.pol is not a struct
    ss = base_struct;
    ss.pol = pi;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol is not a struct', testFailed, verbose);

    % ss.pol.grids is not a struct
    ss = base_struct;
    ss.pol.grids = "grid";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol.grids is not a struct', testFailed, verbose);

    % ss.pol.values is not a struct
    ss = base_struct;
    ss.pol.values = "values";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.pol.values is not a struct', testFailed, verbose);

    % ss.d is not a struct
    ss = base_struct;
    ss.d = 1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.d is not a struct', testFailed, verbose);

    % ss.d.grids is not a struct
    ss = base_struct;
    ss.d.grids = "grid";
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.d.grids is not a struct', testFailed, verbose);

    % ss.agg is not a struct
    ss = base_struct;
    ss.agg = 0;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'ss.agg is not a struct', testFailed, verbose);

    %% === SHOCKS TESTS ===

    % Missing `shocks` field
    ss = base_struct;
    ss = rmfield(ss, 'shocks');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `shocks` field', testFailed, verbose);

    % Missing `shocks.grids`
    ss = base_struct;
    ss.shocks = rmfield(ss.shocks, 'grids');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `shocks.grids`', testFailed, verbose);

    % Wrong type for `shocks.grids.eps_e`
    ss = base_struct;
    ss.shocks.grids.eps_e = ones(3,3);  % invalid: cell array
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'shocks.grids.eps_e are not dense real vectors', testFailed, verbose);

    % Missing both `shocks.Pi` and `shocks.w`
    ss = base_struct;
    ss.shocks = rmfield(ss.shocks, 'w');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'No discretization method (Pi or w) provided', testFailed, verbose);

    % Both `shocks.Pi` and `shocks.w` provided (mutual exclusivity)
    ss = base_struct;
    ss.shocks.Pi = struct('e', ones(size(ss.pol.grids.e,1)));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Both Pi and w should not coexist', testFailed, verbose);

    % AR(1) Tests
    % Shock grid not in var(heterogeneity=) list
    ss = base_struct;
    ss.shocks.grids.fake_shock = ss.shocks.grids.eps_e;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Shock grid symbol not in heterogeneity declaration', testFailed, verbose);

    % Missing Markov matrix for a shock
    ss = base_struct;
    ss.shocks.w = rmfield(ss.shocks.w, 'eps_e');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing w entry for a declared shock', testFailed, verbose);

    % Wrong type for `shocks.grids.w.e`
    ss = base_struct;
    ss.shocks.w.eps_e = speye(numel(ss.shocks.grids.eps_e));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'shocks.w.eps_e are not dense real vectors', testFailed, verbose);

    % Gauss-Hermite weights count
    ss = base_struct;
    ss.shocks.w.eps_e = rand(sizes.N_e, 1);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        ' row count mismatch with grid', testFailed, verbose);

    % Gauss-Hermite weights have negative elements
    ss = base_struct;
    ss.shocks.w.eps_e(1) = -0.1;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Negative elements in Pi', testFailed, verbose);

    % Gauss-Hermite weights do not sum to 1
    ss = base_struct;
    ss.shocks.w.eps_e = ss.shocks.w.eps_e * 2;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Rows of Pi not summing to 1', testFailed, verbose);


    %% === POLICY TESTS ===

    % Missing `pol` field
    ss = base_struct;
    ss = rmfield(ss, 'pol');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol` field', testFailed, verbose);

    % Missing `pol.grids`
    ss = base_struct;
    ss.pol = rmfield(ss.pol, 'grids');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol.grids`', testFailed, verbose);

    % Missing one state grid
    ss = base_struct;
    ss.pol.grids = rmfield(ss.pol.grids, 'a');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing state grid for a', testFailed, verbose);

    % Wrong type for `pol.grids.a`
    ss = base_struct;
    ss.pol.grids.a = ones(3, 3);  % should be a vector
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.grids.a should be a vector', testFailed, verbose);

    % Missing `pol.values`
    ss = base_struct;
    ss.pol = rmfield(ss.pol, 'values');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `pol.values`', testFailed, verbose);

    % Missing policy value for declared variable
    ss = base_struct;
    ss.pol.values = rmfield(ss.pol.values, 'a');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing policy value for a', testFailed, verbose);

    % Wrong type for `pol.values.a`
    ss = base_struct;
    ss.pol.values.a = rand(2, 2, 2); 
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.values.a should be a matrix', testFailed, verbose);

    % Incompatible policy matrix row size (shocks × states)
    ss = base_struct;
    ss.pol.values.a = rand(99, size(ss.pol.values.a, 2));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of rows in policy matrix', testFailed, verbose);

    % Incompatible policy matrix column size
    ss = base_struct;
    ss.pol.values.a = rand(size(ss.pol.values.a, 1), 99);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of columns in policy matrix', testFailed, verbose);

    %% === DISTRIBUTION TESTS ===

    % Missing `d` field
    ss = base_struct;
    ss = rmfield(ss, 'd');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing distribution structure `d`', testFailed, verbose);

    % Wrong type for `ss.d.grids.a`
    ss = base_struct;
    ss.d.grids.a = rand(3, 3);  % not a vector
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.grids.a should be a vector', testFailed, verbose);

    % Missing `d.hist`
    ss = base_struct;
    ss.d = rmfield(ss.d, 'hist');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing d.hist', testFailed, verbose);

    % Wrong type for `ss.d.hist`
    ss = base_struct;
    ss.d.hist = rand(2,2,2);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.hist should be a matrix', testFailed, verbose);

    % Wrong row size in `d.hist`
    ss = base_struct;
    ss.d.hist = rand(99, size(ss.d.hist,2));
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of rows in d.hist', testFailed, verbose);

    % Wrong column size in `d.hist`
    ss = base_struct;
    ss.d.hist = rand(size(ss.d.hist,1), 99);
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Wrong number of columns in d.hist', testFailed, verbose);

    %% === AGGREGATES TESTS ===

    % Missing `agg` field
    ss = base_struct;
    ss = rmfield(ss, 'agg');
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing `agg` field', testFailed, verbose);

    % Remove one required aggregate variable
    ss = base_struct;
    ss.agg = rmfield(ss.agg, 'r');  % remove interest rate
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'Missing aggregate variable r', testFailed, verbose);

    % Wrong type for `ss.agg.r`
    ss = base_struct;
    ss.agg.r = 'non-numeric';
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'agg.r should be numeric', testFailed, verbose);

    %% === PERMUTATION TESTS ===

    % pol.shocks is not a string or cellstr
    ss = base_struct;
    ss.pol.shocks = 123;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.shocks is not a string array or cellstr', testFailed, verbose);

    % pol.shocks has wrong names
    ss = base_struct;
    ss.pol.shocks = {'wrong_name'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.shocks not matching exo_names', testFailed, verbose);

    % pol.states is not a string or cellstr
    ss = base_struct;
    ss.pol.states = 42;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.states is not a string array or cellstr', testFailed, verbose);

    % pol.states has wrong names
    ss = base_struct;
    ss.pol.states = {'bad_state'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'pol.states not matching state_var', testFailed, verbose);

    % d.shocks is not a string or cellstr
    ss = base_struct;
    ss.d.shocks = 3.14;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.shocks is not a string array or cellstr', testFailed, verbose);

    % d.shocks inconsistent with pol.shocks
    ss = base_struct;
    ss.d.shocks = {'ghost_shock'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.shocks not matching pol.shocks', testFailed, verbose);

    % d.states is not a string or cellstr
    ss = base_struct;
    ss.d.states = true;
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.states is not a string array or cellstr', testFailed, verbose);

    % d.states inconsistent with pol.states
    ss = base_struct;
    ss.d.states = {'phantom_state'};
    testFailed = expect_error(@() hank.check_steady_state_input(M_, options_, ss), ...
        'd.states not matching pol.states', testFailed, verbose);

    %% === REDUNDANCY TESTS ===

    options_.hank.nowarningredundant = false;

    % shocks.w contains redundant entry (not in shocks.grids)
    ss = base_struct;
    ss.shocks.xtra = ones(length(ss.shocks.grids.eps_e));  % not listed in shocks.grids
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant shocks.Pi entry accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant shocks.Pi entry: %s\n', ME.message);
    end

    % pol.values contains redundant field (not in model's pol_symbs)
    ss = base_struct;
    ss.pol.values.redundant_var = ones(size(ss.pol.values.a));
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant pol.values field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant pol.values field: %s\n', ME.message);
    end

    % pol.values contains valid optional multiplier field (e.g. MULT_L_a)
    ss = base_struct;
    ss.pol.values.MULT_L_a = 0.1 * ones(size(ss.pol.values.a));
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Optional multiplier MULT_L_a accepted in pol.values (no redundant variable warning expected)');
    catch ME
        testFailed = testFailed + 1;
        fprintf('❌ Unexpected error when supplying pol.values.MULT_L_a: %s\n', ME.message);
    end

    % d.grids contains redundant field (not a declared state)
    ss = base_struct;
    ss.d.grids.not_a_state = [0;1;2];  % not in M_.heterogeneity.state_var
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant d.grids field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant d.grids field: %s\n', ME.message);
    end

    % agg contains extra field (not in M_.endo_names)
    ss = base_struct;
    ss.agg.extra_agg = 999;
    try
        [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
        disp('✔ Redundant agg field accepted (warning expected)');
    catch ME
        testFailed = testFailed+1;
        fprintf('❌ Unexpected error from redundant agg field: %s\n', ME.message);
    end

    if testFailed > 0
        error('Some unit tests associated with the routine `hank.check_steady_state_input` failed!');
    end

end;