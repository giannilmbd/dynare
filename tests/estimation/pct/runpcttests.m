if ~exist('OCTAVE_VERSION', 'builtin') && ~verLessThan('matlab', '23.1')  && matlab.internal.parallel.isPCTInstalled
    f = findobj('Tag', 'dynare-figure');
    delete(f);
    clear;

    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.* ; % in newer MATLAB "import matlab.unittest.plugins.codecoverage.CoverageResult" would be sufficient, but creates lint error in R2020a

    source_dir = getenv('source_root');
    addpath([source_dir filesep 'matlab']);


    suite = testsuite("tParallel");

    runner = testrunner("textoutput");
    format = CoverageResult; %available after R2023a
    % This only covers the main MATLAB process, but that should be enough
    p = CodeCoveragePlugin.forFile(fullfile(source_dir,'matlab','estimation','posterior_sampler_core.m'),'Producing',format); % to prevent lint error in older MATLAB, use old name-value syntax
    runner.addPlugin(p)

    runner.run(suite);

    result = format.Result;
    summary = coverageSummary(result,"statement");
    assert(summary(1)/summary(2) > 0.7, 'Not enough coverage')
end

% generateHTMLReport(result)