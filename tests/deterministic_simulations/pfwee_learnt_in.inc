/* Tests perfect_foresight_with_expectation_errors_{setup,solver}
   using the shocks(learnt_in=…), mshocks(learnt_in=…) and endval(learnt_in=…) syntax
   with dates (instead of integer indices) for periods */

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

check;

// p contains period indices, either as integer or as dates
@#if dates
@# define p = ["2018Y", "2019Y", "2020Y", "2021Y", "2022Y", "2023Y", "2024Y", "2025Y" ]
@#else
@# define p = 1:8
@#endif

@#if bare_first_info_period
shocks;
@#else
shocks(learnt_in = @{p[1]});
@#endif
  var x;
  periods @{p[1]}:@{p[2]};
  values 1.2;
end;

@#if bare_first_info_period
endval;
@#else
endval(learnt_in = @{p[1]});
@#endif
  x = 1.05;
end;

shocks(learnt_in = @{p[2]});
  var x;
  periods @{p[2]};
  add 0.1;
end;

endval(learnt_in = @{p[2]});
  x = 1.1;
end;

shocks(learnt_in = @{p[3]});
  var x;
  periods @{p[3]};
  values 1.4;
end;

mshocks(learnt_in = @{p[3]});
  var x;
  periods @{p[7]};
  values (1.5/1.2); // 1.2 is the terminal steady as anticipated in period 3
end;

endval(learnt_in = @{p[3]});
  x += 0.1;
end;

// Dummy block, that will be overwritten by the next one
shocks(learnt_in = @{p[6]});
  var x;
  periods @{p[6]}:@{p[8]};
  values 10;
end;

shocks(learnt_in = @{p[6]}, overwrite);
  var x;
  periods @{p[6]}:@{p[7]};
  multiply 0.8;
end;

endval(learnt_in = @{p[6]});
  x *= 0.75;
end;

// Save initial steady state (it will be modified by pfwee)
orig_steady_state = oo_.steady_state;
orig_exo_steady_state = oo_.exo_steady_state;

perfect_foresight_with_expectation_errors_setup(
@#if dates
first_simulation_period = 2018Y, last_simulation_period = 2024Y
@#else
periods = 7
@#endif
);

perfect_foresight_with_expectation_errors_solver;
pfwee_simul = oo_.endo_simul;

// Now compute the solution by hand to verify the results
verbatim;

oo_.steady_state = orig_steady_state;
oo_.exo_steady_state = orig_exo_steady_state;

oo_ = make_ex_(M_,options_,oo_);
oo_ = make_y_(M_,options_,oo_);

@#if dates
% Reset {first,last}_simulation_period since we’re going to manipulate options_.periods
options_.simul.first_simulation_period = dates();
options_.simul.last_simulation_period = dates();
options_.periods = 7;
@#endif

% Information arriving in period 1 (temp shock now and tomorrow + permanent shock in future)
oo_.exo_simul(2:3,1) = 1.2;
oo_.exo_simul(4:end,1) = 1.05;
oo_.exo_steady_state = 1.05;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
oo_=perfect_foresight_solver(M_, options_, oo_);

% Information arriving in period 2 (temp shock now + permanent shock in future)
oo_.exo_simul(3,1) = 1.3;
oo_.exo_steady_state = 1.1;
oo_.exo_simul(4:end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 6;
saved_endo = oo_.endo_simul(:, 1);
saved_exo = oo_.exo_simul(1, :);
oo_.endo_simul = oo_.endo_simul(:, 2:end);
oo_.exo_simul = oo_.exo_simul(2:end, :);
oo_=perfect_foresight_solver(M_, options_, oo_);
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

% Information arriving in period 3 (temp shocks + permanent shock in future)
oo_.exo_simul(4,1) = 1.4;
oo_.exo_steady_state = 1.1+0.1;
oo_.exo_simul(8,1) = (1.5/1.2)*oo_.exo_steady_state;
oo_.exo_simul(5:7, 1) = oo_.exo_steady_state;
oo_.exo_simul(9:end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 5;
saved_endo = oo_.endo_simul(:, 1:2);
saved_exo = oo_.exo_simul(1:2, :);
oo_.endo_simul = oo_.endo_simul(:, 3:end);
oo_.exo_simul = oo_.exo_simul(3:end, :);
oo_=perfect_foresight_solver(M_, options_, oo_);
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

% Information arriving in period 6 (temp shocks + permanent shock)
oo_.exo_simul(7,1) = (1.1+0.1)*0.8;
oo_.exo_simul(8,1) = 1.5*0.8;
oo_.exo_steady_state = (1.1+0.1)*0.75;
oo_.exo_simul(9:end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 2;
saved_endo = oo_.endo_simul(:, 1:5);
saved_exo = oo_.exo_simul(1:5, :);
oo_.endo_simul = oo_.endo_simul(:, 6:end);
oo_.exo_simul = oo_.exo_simul(6:end, :);
oo_=perfect_foresight_solver(M_, options_, oo_);
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

% We should have strict equality with first pfwee simulation, because algorithm
% and guess values are exactly the same.
if any(any(pfwee_simul-oo_.endo_simul ~= 0))
    error('Error in perfect_foresight_with_expectation_errors')
end

@#if dates
if ~isequal(Simulated_time_series.dates, dates('2017Y'):dates('2025Y'))
    error('Incorrect Simulated_time_series object')
end
@#endif

end; // verbatim
