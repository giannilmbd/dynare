classdef tParallel < matlab.unittest.TestCase

    methods (TestClassTeardown)
        function tRemoveFigures(~)
            obj = findobj('Tag', 'dynare-figure');
            delete(obj)
            delete(gcp('nocreate'))
        end
    end
    methods (Test)
        function tFmincon(testCase)
            DynareInfoSerial = dynare('SW_2007_fmincon_parallel','nograph','nolog','-DParallel=false');
            out_serial = loadResults("oo_");

            p = gcp('nocreate');
            if isempty( p )
                p = parpool('Processes');
                cobj = onCleanup(@() delete(p));
            elseif ~isempty(gcp('nocreate')) && isa(p, 'parallel.ThreadPool')
                delete(p)
                p = parpool('Processes');
                cobj = onCleanup(@() delete(p));
            end
            DynareInfoParallel = dynare('SW_2007_fmincon_parallel','nograph','nolog','-DParallel=true');
            out_parallel = loadResults("oo_");

            testCase.verifyEqual(out_serial.oo_.posterior.optimization.log_density, out_parallel.oo_.posterior.optimization.log_density, RelTol = 2e-9)
            testCase.verifyEqual(out_serial.oo_.MarginalDensity.LaplaceApproximation,out_parallel.oo_.MarginalDensity.LaplaceApproximation, RelTol = 2e-7)

            fprintf('Parallel optim with %d workers of %d threads was %0.2f%% faster\n', p.NumWorkers,  p.NumThreads, (1-DynareInfoParallel.time.compute./DynareInfoSerial.time.compute)*100)
        end

        function tParallelFs2000(testCase)

            c = parcluster('Processes');
            c.NumThreads = floor(c.NumWorkers/2);
            p = gcp('nocreate');
            if ~isempty(p)
                delete(p)
            end
            p = parpool(c, 2);
            DynareInfoParallel = dynare('fs2000.mod', 'nograph','nolog');
            out_parallel = loadResults("oo_");

            nWork = p.NumWorkers;
            nThrea = p.NumThreads;
            delete(p);

            DynareInfoSerial = dynare('fs2000.mod', 'nograph','nolog');
            out_serial = loadResults("oo_");

            testCase.verifyEqual(out_serial.oo_.posterior.optimization.log_density, out_parallel.oo_.posterior.optimization.log_density, RelTol = 2e-9)
            testCase.verifyEqual(out_serial.oo_.MarginalDensity.LaplaceApproximation,out_parallel.oo_.MarginalDensity.LaplaceApproximation, RelTol = 2e-7)

            % remove timing fields, these won't match
            out.oo_ = rmfield(out_parallel.oo_, 'time');
            out2.oo_ = rmfield(out_serial.oo_, 'time');

            % Start comparison
            flds1 = fields(out);
            flds2 = fields(out2);
            testCase.verifyEqual(flds1, flds2)


            for i = 1 : numel(flds1)
                testCase.verifyEqual(out.(flds1{i}), out2.(flds1{i}))
            end

            fprintf('Parallel optim with %d workers of %d threads was %0.2f%% faster\n', nWork, nThrea, (1-DynareInfoParallel.time.compute./DynareInfoSerial.time.compute)*100)

        end
    end

end

function out = loadResults(vars)

arguments
    vars (1,:) string {mustBeNonempty}
end

[~,tname] = fileparts(tempname);
tname = tname+".mat";

cmd = sprintf("save %s %s", tname, strjoin(vars, " "));
evalin('base',cmd);
out = load(tname);
delete(tname)

end