/*
 * Contains the core estimation logic used by the serial and PCT test files.
 * It runs estimations with different seeding strategies (start/default, specific seed 123,
 * no seed, reset to default) and performs checks to verify that seeds behave as expected
 * (e.g., specific seeds produce identical results, different seeds produce different results).
 */

/*
 * Copyright © 2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

@#define MH_REPLIC  = 3000

% first run of a mod file uses set_dynare_seed('default');
fprintf('Seed at start:\n');
disp(options_.DynareRandomStreams);
estimation(datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
           mode_compute=0,nograph,nodisplay,
           mode_file='ls2003/Output/ls2003_mode',
           mh_nblocks=@{MH_NBLOCKS},mh_jscale=0.5,mh_replic=@{MH_REPLIC},
           dirname='@{DIRNAME}_start');
@#for i in 1:MH_NBLOCKS
seed_start_@{i} = load(['@{DIRNAME}_start/metropolis/' M_.fname '_mh1_blck@{i}.mat']);
@#endfor

set_dynare_seed(123); % should give different results to default
fprintf('Seed after set_dynare_seed(''123''):\n');
disp(options_.DynareRandomStreams);
estimation(datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
           mode_compute=0,nograph,nodisplay,
           mode_file='ls2003/Output/ls2003_mode',
           mh_nblocks=@{MH_NBLOCKS},mh_jscale=0.5,mh_replic=@{MH_REPLIC},
           dirname='@{DIRNAME}_123');
@#for i in 1:MH_NBLOCKS
seed_123_@{i} = load(['@{DIRNAME}_123/metropolis/' M_.fname '_mh1_blck@{i}.mat']);
@#endfor

fprintf('Seed without any set_dynare_seed command:\n');
disp(options_.DynareRandomStreams);
estimation(datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
           mode_compute=0,nograph,nodisplay,
           mode_file='ls2003/Output/ls2003_mode',
           mh_nblocks=@{MH_NBLOCKS},mh_jscale=0.5,mh_replic=@{MH_REPLIC},
           dirname='@{DIRNAME}_noseed');
@#for i in 1:MH_NBLOCKS
noseed_@{i} = load(['@{DIRNAME}_noseed/metropolis/' M_.fname '_mh1_blck@{i}.mat']);
@#endfor

set_dynare_seed('default'); % should be equal to serial_start
fprintf('Seed after set_dynare_seed(''default''):\n');
disp(options_.DynareRandomStreams);
estimation(datafile=data_ca1,first_obs=8,nobs=79,prefilter=1,
           mode_compute=0,nograph,nodisplay,
           mode_file='ls2003/Output/ls2003_mode',
           mh_nblocks=@{MH_NBLOCKS},mh_jscale=0.5,mh_replic=@{MH_REPLIC},
           dirname='@{DIRNAME}_default');
@#for i in 1:MH_NBLOCKS
seed_default_@{i} = load(['@{DIRNAME}_default/metropolis/' M_.fname '_mh1_blck@{i}.mat']);
@#endfor

@#for i in 1:MH_NBLOCKS

fprintf('\nRESULTS FOR BLOCK @{i}:\n');

msg = sprintf('seed at start and set_dynare_seed(''default'') should be equal!');
if isequaln(seed_start_@{i}.x2, seed_default_@{i}.x2)==1 && isequaln(seed_start_@{i}.logpo2, seed_default_@{i}.logpo2)==1
    fprintf('✅ %s\n',msg);
else
    error('❌ %s',msg);
end

msg = sprintf('seed at start and set_dynare_seed(123) should not be equal!');
if isequaln(seed_start_@{i}.x2, seed_123_@{i}.x2)==1 && isequaln(seed_start_@{i}.logpo2, seed_123_@{i}.logpo2)==1
    error('❌ %s',msg);
else
    fprintf('✅ %s\n',msg);
end

msg = sprintf('seed at start and noseed should not be equal!');
if isequaln(seed_start_@{i}.x2, noseed_@{i}.x2)==1 && isequaln(seed_start_@{i}.logpo2, noseed_@{i}.logpo2)==1
    error('❌ %s',msg);
else
    fprintf('✅ %s\n',msg);
end

msg = sprintf('set_dynare_seed(''default'') and set_dynare_seed(123) should not be equal!');
if isequaln(seed_default_@{i}.x2, seed_123_@{i}.x2)==1 && isequaln(seed_default_@{i}.logpo2, seed_123_@{i}.logpo2)==1
    error('❌ %s',msg);
else
    fprintf('✅ %s\n',msg);
end

msg = sprintf('set_dynare_seed(''default'') and noseed (which came directly after set_dynare_seed(''123'') should not be equal!');
if isequaln(seed_default_@{i}.x2, noseed_@{i}.x2)==1 && isequaln(seed_default_@{i}.logpo2, noseed_@{i}.logpo2)==1
    error('❌ %s',msg);
else
    fprintf('✅ %s\n',msg);
end

msg = sprintf('set_dynare_seed(123) and noseed (which came directly after set_dynare_seed(''123'')) should be equal!');
if isequaln(seed_123_@{i}.x2, noseed_@{i}.x2)==1 && isequaln(seed_123_@{i}.logpo2, noseed_@{i}.logpo2)==1
    fprintf('✅ %s\n',msg);
else
    error('❌ %s',msg);
end

@#endfor
fprintf('\n');
